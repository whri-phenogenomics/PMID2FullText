import logging
from dataclasses import dataclass
from pathlib import Path

from Bio import Entrez
import click

from src.PMID2FullText.agents.excel_agent import ExcelAgent

Entrez.email = 'c.kroll@qmul.ac.uk'

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


@dataclass
class NCBIAgent:
    no_pmcid_output_file: str

    def __post_init__(self):
        self.no_pmcid_output_file = self.resolve_path(self.no_pmcid_output_file)

    def resolve_path(self, path):
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
            self.handle_no_pmcid(pmid)

    def handle_no_pmcid(self, pmid):
        """Write every pmid that is not in pmcid to outfile"""
        with open(self.no_pmcid_output_file, "a") as f:
            f.write(pmid + "\n")

    def is_open_access(self, pmcid):
        """Check if the article is open access."""
        handle = Entrez.esummary(db="pmc", id=pmcid)
        summary = Entrez.read(handle)
        handle.close()
        oa_status = summary[0].get('OpenAccess')
        return oa_status == '1'

    # def download_full_text(self, pmcid, pmid):
    #     """Download the full text of the article."""
    #     handle = Entrez.efetch(db="pmc", id=pmcid, rettype="full", retmode="xml")
    #     full_text = handle.read()
    #     handle.close()
    #     with open(f'{pmid}.xml', 'w') as f:
    #         f.write(full_text)
    #
    # def process_pmid(self, pmid):
    #     """Process each PMID: check availability and download if possible."""
    #     pmcid = self.get_pmcid(pmid)
    #     if pmcid:
    #         if self.is_open_access(pmcid):
    #             self.download_full_text(pmcid, pmid)
    #             print(f'PMID {pmid}: Full text downloaded.')
    #         else:
    #             print(f'PMID {pmid}: Full text not open access.')
    #     else:
    #         print(f'PMID {pmid}: No PMCID found, full text not available.')

    def check(self, pmid: str) -> bool:
        """Check if the PMID has a full text available (open access)."""
        pmcid = self.get_pmcid(pmid)
        if pmcid:
            return self.is_open_access(pmcid)
        else:
            self.handle_no_pmcid(pmid)
            return False


@click.command()
@click.option('-f', '--file_path', default='input/data.xlsx', help="Path to the file containing PMIDs.")
@click.option('-o', '--output-dir', default="output/no_pmcids.txt" ,help="Directory to save XML files.")
@click.option('--verbose', is_flag=True, help="Enable verbose logging.")
def main(file_path, output_dir, verbose):
    if verbose:
        logger.setLevel(logging.INFO)
    excel = ExcelAgent(path=file_path)
    pmids = excel.read_pmids_from_xlsx(column_index=9)
    ncbi_agent = NCBIAgent(no_pmcid_output_file=output_dir)
    for pmid in pmids:
        ncbi_agent.check(pmid)
    logger.info("Process completed successfully.")

if __name__ == '__main__':
    main()
