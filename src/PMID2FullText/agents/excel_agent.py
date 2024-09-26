import logging
from dataclasses import dataclass
from pathlib import Path
from typing import List

import pandas as pd

logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


@dataclass
class ExcelAgent:
    path: str

    def __post_init__(self):
        self.path = self.resolve_path(self.path)
        logger.info(self.path)

    def resolve_path(self, path):
        project_root = Path(__file__).parent.parent
        full_path = project_root / path
        if not full_path.exists():
            raise FileNotFoundError(f"File not found: {full_path}")
        return str(full_path)

    def read_pmids_from_xlsx(self, column_index=9) -> List[str]:
        df = pd.read_excel(self.path, dtype=str)  # read as str to avoid .o at end of numbers
        return df.iloc[:, column_index].dropna().str.strip().tolist()
