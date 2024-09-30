import os
from pathlib import Path

from bs4 import BeautifulSoup
from dataclasses import dataclass
from typing import List, Tuple


@dataclass
class PubmedXMLWrapper:
    directory: str

    def __post_init__(self):
        self.directory = self.resolve_path(self.directory)
        print(self.directory)

    @staticmethod
    def resolve_path(path):
        project_root = Path(__file__).parent.parent
        full_path = project_root / path
        return str(full_path)

    def parse_pmxml(self, xml_content: str) -> BeautifulSoup:
        """Parse the given XML content into a BeautifulSoup object."""
        return BeautifulSoup(xml_content, "xml")

    def check_body_tag(self, soup: BeautifulSoup) -> bool:
        """Check if the given BeautifulSoup object has a <body> tag."""
        try:
            return bool(soup.find("pmc-articleset").find("article").find("body"))
        except AttributeError:
            return False

    def get_title(self, soup: BeautifulSoup) -> str:
        """Retrieve the article title from the given BeautifulSoup object."""
        title_tag = soup.find("ArticleTitle")
        return title_tag.text if title_tag else "No title available"

    def get_abstract(self, soup: BeautifulSoup) -> str:
        """Retrieve the abstract from the given BeautifulSoup object."""
        abstract_tag = soup.find("Abstract")
        return abstract_tag.text if abstract_tag else "No abstract available"

    def process_files(self) -> Tuple[int, int, List[str]]:
        """Process all XML files in the directory and count how many have no <body> tag."""
        no_body_count = 0
        total_files = 0
        no_body_files = []

        for file_name in os.listdir(self.directory):
            if file_name.endswith(".xml"):
                total_files += 1
                file_path = os.path.join(self.directory, file_name)

                with open(file_path, 'r', encoding='utf-8') as xml_file:
                    xml_content = xml_file.read()

                soup = self.parse_pmxml(xml_content)

                # Check if the <body> tag exists
                if not self.check_body_tag(soup):
                    no_body_count += 1
                    no_body_files.append(file_name)

                # title = self.get_title(soup)
                # abstract = self.get_abstract(soup)
                # print(f"File: {file_name}\nTitle: {title}\nAbstract: {abstract}\n")

        return total_files, no_body_count, no_body_files


if __name__ == "__main__":
    path = Path("output_final")
    xml_wrapper = PubmedXMLWrapper(directory=str(path))
    total, no_body, no_body_files = xml_wrapper.process_files()

    print(f"Total files processed: {total}")
    print(f"Files without <body> tag: {no_body}")
    if no_body_files:
        print("Files missing <body> tag:")
        for file_name in no_body_files:
            print(file_name)
