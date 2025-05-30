from pathlib import Path
import pandas as pd

path = Path(__file__).parent

standards_path = path / "standards.txt"
standards_df = pd.read_csv(
    standards_path,
    names=['name', 'ra', 'dec', 'format', 'type'],
    sep='\s\s+', engine='python'
)
