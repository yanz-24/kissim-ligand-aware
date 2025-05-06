"""
kissim.io.dataframe

Defines a DataFrame-based pocket class.
"""
import pandas as pd

from opencadd.structure.pocket import Pocket

from . import KlifsToKissimData


class PocketDataFrame(Pocket):
    @classmethod
    def from_structure_klifs_id(cls, structure_klifs_id, klifs_session=None):
        """
        Get DataFrame-based pocket object from a KLIFS structure ID.

        Parameters
        ----------
        structure_id : int
            KLIFS structure ID.
        klifs_session : None or opencadd.databases.klifs.session.Session
            Local or remote KLIFS session. If None, generate new remote session.

        Returns
        -------
        kissim.io.PocketDataFrame or None
            DataFrame-based pocket object.
        """

        data = KlifsToKissimData.from_structure_klifs_id(structure_klifs_id, klifs_session)
        if data:
            pocket = cls.from_text(
                data.text, data.extension, data.residue_ids, data.residue_ixs, structure_klifs_id, data.ligand_expo_id,
            )
            return pocket
        else:
            return None
        

class LigandDataFrame():
    """
    DataFrame-based object for ligand data.
    """
    def __init__(self, ligand_expo_id: str):
        """
        Initialize the LigandDataFrame object.
        Parameters
        ----------
        ligand_expo_id : str
            Ligand expo ID.
        ligand_df : pd.DataFrame
            DataFrame containing ligand data.
        """
        self.df = pd.DataFrame()
        self.ligand_expo_id = ligand_expo_id.upper()
         
    @classmethod
    def from_text(cls, text: str, extension: str, ligand_expo_id: str):
        """
        Extract HETATM information from PDB text for specified ligand_expo_id
        """
        if extension != "pdb":
            raise ValueError(f"Unsupported file extension: {extension}. Only 'pdb' is supported.")
        lines = text.splitlines()
        records = []
        ligand_expo_id = ligand_expo_id.upper()

        for line in lines:
            if line.startswith("HETATM") and line[17:20].strip().upper() == ligand_expo_id:
                atom_id = int(line[6:11])
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip()
                residue_id = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                records.append({
                    'atom.id': atom_id,
                    'atom.name': atom_name,
                    'atom.x': x, 'atom.y': y, 'atom.z': z,
                    'residue.id': residue_id,
                    'residue.name': residue_name,
                })

        if not records:
            raise ValueError(f"No atoms found for ligand_expo_id='{ligand_expo_id}'")

        df = pd.DataFrame(records)
        obj = cls(ligand_expo_id)
        obj.df = df
        return obj


    @classmethod
    def from_structure_klifs_id(cls, structure_klifs_id, klifs_session=None):
        """
        Get DataFrame-based ligand object from a KLIFS ligand ID.

        Parameters
        ----------
        ligand_expo_id : str
            KLIFS ligand ID.
        klifs_session : None or opencadd.databases.klifs.session.Session
            Local or remote KLIFS session. If None, generate new remote session.

        Returns
        -------
        kissim.io.LigandDataFrame or None
            DataFrame-based ligand object.
        """
        data = KlifsToKissimData.from_structure_klifs_id(structure_klifs_id, klifs_session)
        if data:
            ligand = cls.from_text(
                data.text, data.extension, data.ligand_expo_id
            )
            return ligand
        else:
            return None
            
    
