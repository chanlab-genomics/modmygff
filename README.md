# modmygff

## Usage

`modmygff` is a tool to add (possibly multiple) external database references to an existing gff file. It's usage is best shown with a simple example. Take here this gff file (`Polarella_glacialis_CCMP2088.gff3`) with the following layout
```
CCMP2088_scf7180000719957 | EVM	gene	13904	57986	.	-	.	ID=evm.TU.CCMP2088_scf7180000719957.1;Name=evm.model.CCMP2088_scf7180000719957.1

CCMP2088_scf7180000719957	EVM	mRNA	13904	57986	.	-	.	ID=evm.model.CCMP2088_scf7180000719957.1;Parent=evm.TU.CCMP2088_scf7180000719957.1;Name=evm.model.CCMP2088_scf7180000719957.1

CCMP2088_scf7180000719957	EVM	exon	57958	57986	.	-	.	ID=evm.model.CCMP2088_scf7180000719957.1.exon1;Parent=evm.model.CCMP2088_scf7180000719957.1

CCMP2088_scf7180000719957	EVM	CDS	57958	57986	.	-	0	ID=cds.evm.model.CCMP2088_scf7180000719957.1;Parent=evm.model.CCMP2088_scf7180000719957.1

CCMP2088_scf7180000719957	EVM	exon	54190	54247	.	-	.	ID=evm.model.CCMP2088_scf7180000719957.1.exon2;Parent=evm.model.CCMP2088_scf7180000719957.1
```
Suppose we would like to add the following add information from two tabular separated containing external database references. The first external database reference dataframe is held in `CCMP2088_UniProt.tsv` and has the following layout

|                                       |                                  |
| ------------------------------------- | -------------------------------- |
| evm.model.CCMP2088_scf7180000666762.1 | tr\|G1KF44\|G1KF44_ANOCA         |
| evm.model.CCMP2088_scf7180000666762.2 | tr\|A0A0G4FFG2\|A0A0G4FFG2_VITBC |
| evm.model.CCMP2088_scf7180000666765.1 | tr\|A0A1Q9E259\|A0A1Q9E259_SYMMI |
| evm.model.CCMP2088_scf7180000666765.2 | tr\|A0A1Q9ESA1\|A0A1Q9ESA1_SYMMI |
| evm.model.CCMP2088_scf7180000666766.1 | tr\|C5LC67\|C5LC67_PERM5         |

Our other dataframe is held in `CCMP2088_pfam.tsv` as has the following layout

|                                       |     |     |     |     |            |
| ------------------------------------- | --- | --- | --- | --- | ---------- |
| evm.model.CCMP2088_scf7180000720760.2 | 412 | 688 | 408 | 690 | PF07714.16 |
| evm.model.CCMP2088_scf7180000720760.3 | 196 | 604 | 191 | 604 | PF00501.27 |
| evm.model.CCMP2088_scf7180000720760.3 | 612 | 677 | 612 | 677 | PF13193.5  |
| evm.model.CCMP2088_scf7180000720760.4 | 281 | 314 | 274 | 320 | PF07973.13 |
| evm.model.CCMP2088_scf7180000720760.5 | 13  | 85  | 12  | 85  | PF01423.21 |

Now is a good time to mention that this program will only accept external database reference dataframe that are *tabular separated*. Notice that in `CCMP2088_UniProt.tsv` the gene ID is held in the zeroth column (zero indexed) while the references are held in the first column and for `CCMP2088_pfam.tsv` gene ID is again held in the zeroth column (zero indexed) while the references are instead held in the fifth column. We now have everything we need to add reference from both these files to our gff3 file. The `modmygff.py` script (which is used to drive the modifying) will expect three command line arguments:
- `--gff_path` This is just a file path to the gff file that will have the database references added to it.
- `--annotation` The annotation argument will expect three arguments to come after it. The first being the path to the annotation dataframe which is holding the database references. The other two argument in the (numerical, zero-indexed) position of the gene ID column and database reference column within the tabular separated file.
- `--output_path` This is the output path of the modified gff file.

To run the above example we could use
```
python .\modmygff.py --gff_path "Polarella_glacialis_CCMP2088.gff3" --annotation "CCMP2088_UniProt.tsv" 0 1  --annotation "CCMP2088_pfam.tsv" 0 5 --output_path Polarella_glacialis_CCMP2088_ext.gff3"
```
Once running `modmygff.py` should show you its progress
```
Reading gff file
ERROR    Line 1: FORMAT: "##gff-version" missing from the first line
-> CCMP2088_scf7180000719957    EVM     gene    13904   57986   .       -       .       ID=evm.TU.CCMP2088_scf7180000719957.1;Name=evm.model.CCMP2088_scf7180000719957.1
Modify Compilation:  16%|#####4                            | 197225/1224532 [00:16<01:28, 11628.37it/s]
```