"""Pubmed Client."""

import logging
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import List, Union
from urllib import parse
from dotenv import load_dotenv

import click
import requests
from bs4 import BeautifulSoup
from requests import Request

from src.PMID2FullText.agents.excel_agent import ExcelAgent
from src.PMID2FullText.agents.ncbi_agent import NCBIAgent

PMID = str
TITLE_WEIGHT = 5
MAX_PMIDS = 50
RETRY_MAX = 5

PUBMED = "pubmed"
EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"



@dataclass
class PubmedClient:
    """A client for the Pubmed API.

    This class is a wrapper around the Entrez API.
    """
    load_dotenv()
    max_text_length: int = 100000000000

    try:
        email = os.environ.get("NCBI_EMAIL")
        logging.info(f"NCBI MAIL found: {email}")
    except ValueError:
        email = None
        logging.info("Email for NCBI API not found.")

    try:
        ncbi_key = os.environ.get("NCBI_API_KEY")
        logging.info(f"NCBI API KEY found: {ncbi_key}")
    except ValueError:
        ncbi_key = None
        logging.info("NCBI API key not found. Will use no key.")


    def text(
        self, ids: Union[list[PMID], PMID], output, raw=False, autoformat=True, pubmedcental=True
    ):
        """Get the text of one or more papers from their PMIDs.

        :param ids: List of PubMed IDs, or string with single PMID
        :param raw: if True, do not parse the xml, just return the raw output with tags
        :param autoformat: if True include title and abstract concatenated
        :param pubmedcentral: if True, retreive text from PubMed Central where possible
        :return: the text of a single entry, or a list of strings for text of multiple entries
        """
        xml_datas = []
        for id_ in ids:
            xml_data = self.pmc_text(id_)
            xml_datas.append(xml_data)
        output_dir = Path(output)
        output_dir.mkdir(parents=True, exist_ok=True)
        for id, xml_text in zip(ids, xml_datas):
            output_file = f"{output_dir}" + f"{id}.xml"
            logging.info(f"Output File Path: {output_file}")
            with open(output_file, 'w', encoding='utf-8') as f:
                f.write(xml_text)
            logging.info(f'Saved XML for PMID {id} to {output_file}')


    def pmc_text(self, pmc_id: str) -> str:
        """Get the text of one PubMed Central entry.

        Don't parse further here - just get the raw response.
        :param pmc_id: List of PubMed IDs, or string with single PMID
        :return: the text of a single entry as XML
        """
        xml_data = ""

        fetch_url = EUTILS_URL + "efetch.fcgi"
        logging.info(f"FETCH URL: {fetch_url}")

        if self.email and self.ncbi_key:
            params = {
                "db": "pmc",
                "id": pmc_id,
                "rettype": "xml",
                "retmode": "xml",
                "email": self.email,
                "api_key": self.ncbi_key,
            }
        else:
            params = {"db": "pmc", "id": pmc_id, "rettype": "xml", "retmode": "xml"}

        req = Request('GET', fetch_url, params=params)
        prepared = req.prepare()
        full_url = prepared.url
        logging.info(f"Requesting URL: {full_url}")  # This logs the full URL

        response = requests.get(fetch_url, params=parse.urlencode(params, safe=","))
        logging.info(f"Response: {response}")

        trying = True
        try_count = 0
        while trying:
            if response.status_code == 200:
                xml_data = response.text
                trying = False
            else:
                logging.error(
                    f"Encountered error in fetching from PubMed Central: {response.status_code}"
                )
                try_count = try_count + 1
                if try_count < RETRY_MAX:
                    logging.info("Trying again...")
                    time.sleep(2)
                else:
                    logging.info(f"Giving up - last status code {response.status_code}")
                    trying = False
        logging.info(f"Retrieved PubMed Central document data for {pmc_id}.")

        return xml_data


    def parse_pmxml(self, xml: str, raw: bool, autoformat: bool, pubmedcentral: bool):
        """Extract structured text from PubMed and PubMed Central XML.

        :param xml: One or more xml entries, as string
        :param raw: if True, do not parse the xml beyond separating documents
        :param autoformat: if True include title and abstract concatenated
            Otherwise the output will include ALL text contents besides XML tags
        :param pubmedcentral: if True replace abstract with PubMed Central text
            If there isn't a PMC ID, just use the abstract.
            If there is a PMC ID, use the abstract AND chunk the body text.
            This means the same ID may have multiple entries and may require multiple
            queries to the LLM.
        :return: a list of strings, one per entry
        """
        xml_data = ""
        soup = BeautifulSoup(xml, "xml")

        logging.info("Parsing all xml entries...")
        for pa in soup.find_all(["PubmedArticle", "PubmedBookArticle"]):
            # First check the PMID, and if requested, any PMC ID
            pmid = ""
            if pa.find("PMID"):
                pmid = pa.find("PMID").text
            pmc_id = ""
            try:  # There's a chance that this entry is missing one or more fields below
                if (
                    pa.find("PubmedData").find("ArticleIdList").find("ArticleId", {"IdType": "pmc"})
                    and pubmedcentral
                ):
                    pmc_id = (
                        pa.find("PubmedData")
                        .find("ArticleIdList")
                        .find("ArticleId", {"IdType": "pmc"})
                        .text
                    )
            except AttributeError:
                logging.info(f"PubMed entry {pmid} is missing the expected PubMedData fields.")

            xml_data = self.pmc_text(pmc_id)
        return xml_data


@click.command()
@click.option('-f', '--file_path', default= "input/data.xlsx", required=True, help="Path to the file containing PMIDs.")
@click.option('-o', '--output_dir', default='output_fun/pmc_xml', help="Directory to save XML files.")
@click.option('--max-text-length', default=0, type=int, help="Maximum length of text to include in each output file.")
@click.option('--no-pmcid', default="output/no_pmcid.txt", help="File to save PMIDs without PMCID.")
@click.option('--raw', is_flag=True, help="If set, do not parse the XML, just return the raw output with tags.")
@click.option('--verbose', is_flag=True, help="Enable verbose logging.")
def main(file_path, output_dir, max_text_length, no_pmcid, raw, verbose):
    """Download PubMed articles as XML files."""

    if verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)
    agent = NCBIAgent(xml_output_dir=output_dir, no_pmcid_output_file=no_pmcid )
    excel = ExcelAgent(path=file_path)
    pmids = excel.read_pmids_from_xlsx(column_index=9)
    pmcids = []
    for id_ in pmids:
        print(id_)
        pmcid = agent.get_pmcid(id_)
        print(pmcid)
        pmcids.append(pmcid)
    client = PubmedClient(max_text_length=max_text_length)
    client.text(pmcids, output=output_dir)


if __name__ == '__main__':
    main()