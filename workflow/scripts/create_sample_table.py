from pathlib import Path
import pandas as pd

data = []

with open('../config/path_list.txt', 'r') as fp:
    lines = fp.readlines()

for line in lines:
    p = Path(line.strip())
    data.append((p.stem, str(p)))


samples = (
    pd.DataFrame.from_records(data, columns=['sample_name', 'sample_path']).set_index("sample_name", drop=False).sort_index()
)



# Writing to file
# with open("myfile.txt", "w") as fp:
#     fp.writelines(samples)
