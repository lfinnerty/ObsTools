"""
Fetch RA, Dec, and K-band magnitude for a list of TESS Input Catalog (TIC)
IDs using astroquery's MAST Catalogs interface (TIC-8 / Stassun et al. 2019).

Candidates are from:
Cullen & Bayliss (2024), MNRAS 531, 1133 - "A search for non-transiting
exoplanets with optical light phase curves from TESS Southern ecliptic
sectors" (Table 2/3 lists TIC IDs; RA/Dec/Kmag are not given in the paper
itself and are pulled here from the TIC-8 catalog).

Requirements:
    pip install astroquery astropy pandas
"""

from astroquery.mast import Catalogs
import pandas as pd

# TIC IDs for the 27 non-transiting hot Jupiter candidates
tic_ids = [
    124280718, 141372241, 266784171, 351601347, 243494729,
    100512121, 362086194, 62078858, 96918158, 251855019,
    200526405, 2758451, 196322336, 264903281, 174001896,
    258914469, 121026156, 380914081, 56126064, 235055610,
    454198279, 59534077, 144305370, 52199183, 388496589,
    158716775, 443857085,
]

def fetch_tic_data(tic_id):
    """Query MAST TIC-8 catalog for a single TIC ID and return key fields."""
    try:
        result = Catalogs.query_criteria(catalog="Tic", ID=tic_id)
        if len(result) == 0:
            return {"TIC_ID": tic_id, "RA": None, "Dec": None, "Kmag": None,
                    "error": "No match found"}
        row = result[0]
        return {
            "TIC_ID": tic_id,
            "RA": row["ra"],          # degrees, J2000
            "Dec": row["dec"],        # degrees, J2000
            "Kmag": row["Kmag"],      # 2MASS K-band magnitude
            "Tmag": row["Tmag"],      # TESS magnitude, for cross-check
            "error": None,
        }
    except Exception as e:
        return {"TIC_ID": tic_id, "RA": None, "Dec": None, "Kmag": None,
                "error": str(e)}


def main():
    rows = [fetch_tic_data(tid) for tid in tic_ids]
    df = pd.DataFrame(rows)

    # Report any failures
    failed = df[df["error"].notna()]
    if not failed.empty:
        print("Warning: some lookups failed:")
        print(failed[["TIC_ID", "error"]])

    # Save results
    out_path = "tic_candidate_data.csv"
    df.to_csv(out_path, index=False)
    print(f"\nSaved {len(df)} rows to {out_path}\n")
    print(df.to_string(index=False))


if __name__ == "__main__":
    main()
