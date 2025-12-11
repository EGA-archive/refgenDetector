import pandas as pd
import glob

# Process all pickle files
for pkl_file in glob.glob("./pkls/*.pkl"):
    print(f"Processing {pkl_file}...")
    df = pd.read_pickle(pkl_file)
    
    # Set index if not already set
    if df.index.name != 'position' and 'position' in df.columns:
        df = df.set_index('position')
        df.to_pickle(pkl_file)
        print(f"  ✓ Re-saved with position index")
    else:
        print(f"  - Already indexed or no position column")