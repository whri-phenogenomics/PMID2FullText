import logging
from dataclasses import dataclass
from pathlib import Path

from Bio import Entrez
import click

from src.PMID2FullText.agents.excel_agent import ExcelAgent

Entrez.email = 'c.kroll@qmul.ac.uk'

logger = logging.getLogger(__name__)
console_handler = logging.StreamHandler()
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

@dataclass
class NCBIAgent:
    no_pmcid_output_file: str
    xml_output_dir: str

    def __post_init__(self):
        self.no_pmcid_output_file = self.resolve_path(self.no_pmcid_output_file)
        self.xml_output_dir = self.resolve_path(self.xml_output_dir)
        xml_output_dir_path = Path(self.xml_output_dir) / "pmc_xml"
        xml_output_dir_path.mkdir(parents=True, exist_ok=True)

    @staticmethod
    def resolve_path(path):
        project_root = Path(__file__).parent.parent
        full_path = project_root / path
        return str(full_path)

    def get_pmcid(self, pmid):
        """Retrieve the PMCID for a given PMID."""
        handle = Entrez.elink(dbfrom="pubmed", db="pmc", id=pmid)
        record = Entrez.read(handle)
        handle.close()
        try:
            pmcid = record[0]['LinkSetDb'][0]['Link'][0]['Id']
            return pmcid
        except (IndexError, KeyError):
            logger.info(f"No PMCID for {pmid}")

    def handle_no_pmcid(self, pmid):
        """Write every pmid that is not in pmcid to outfile"""
        with open(self.no_pmcid_output_file, "a") as f:
            f.write(pmid + "\n")

    @staticmethod
    def is_open_access(pmcid):
        """Check if the article is open access."""
        handle = Entrez.esummary(db="pmc", id=pmcid)
        summary = Entrez.read(handle)
        handle.close()
        oa_status = summary[0].get('OpenAccess')
        return oa_status == '1'

    def download_full_text(self, pmcid, pmid):
        """Download the full text of the article."""
        try:
            handle = Entrez.efetch(db="pmc", id=pmcid, rettype="full", retmode="xml")
            full_text = handle.read()
            handle.close()
            if "<body>" not in full_text:
                logger.warning(f"Full text for PMCID {pmcid} (PMID {pmid}) is incomplete. Skipping save.")
                return False
            output_file_path = Path(self.xml_output_dir) / "pmc_xml" / f'{pmid}.xml'
            with open(output_file_path, 'w') as f:
                f.write(full_text)
            logger.info(f"Full text for PMID {pmid} (PMCID {pmcid}) downloaded successfully.")
            return True
        except Exception as e:
            logger.error(f"Failed to download full text for PMID {pmid}: {str(e)}")


    def process_pmid(self, pmid):
            """Process each PMID: check availability and download if possible."""
            pmcid = self.get_pmcid(pmid)
            if pmcid:
                if self.is_open_access(pmcid):
                    self.download_full_text(pmcid, pmid)
                    logger.info(f'PMID {pmid}: Full text downloaded.')
                else:
                    logger.info(f'PMID {pmid}: Full text not open access.')
            else:
                logger.info(f'PMID {pmid}: No PMCID found, full text not available.')

    # only to check how many are not on pmc
    def check(self, pmid: str) -> bool:
        """Check if the PMID has a full text available (open access)."""
        pmcid = self.get_pmcid(pmid)
        if pmcid:
            logger.info(f"{self.is_open_access(pmcid)}: {pmid} not open access")
            if self.is_open_access(pmcid):
                logger.info(f'PMID {pmid}: Full text downloadable.')

        else:
            self.handle_no_pmcid(pmid)
            logger.info(f'PMID {pmid}: Full text not open for download.')

            return False


@click.command()
@click.option('-f', '--file_path', default='input/data.xlsx', help="Path to the file containing PMIDs.")
@click.option('-o', '--full_text_xml_output', default="output" ,help="Directory to save XML files.")
@click.option('--no-pmcid', default="output/no_pmcid.txt", help="File to save PMIDs without PMCID.")
@click.option('--verbose', is_flag=True, help="Enable verbose logging.")
def main(file_path, full_text_xml_output, no_pmcid, verbose):
    if verbose:
        logger.setLevel(logging.INFO)
    excel = ExcelAgent(path=file_path)
    pmids = excel.read_pmids_from_xlsx(column_index=9)
    ncbi_agent = NCBIAgent(no_pmcid_output_file=no_pmcid, xml_output_dir=full_text_xml_output)
    for pmid in pmids:
        print("go")
        ncbi_agent.check(pmid)
    logger.info("Process completed successfully.")

if __name__ == '__main__':
    main()
