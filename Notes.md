# Notes

2020/03/06 @ 2:16 PM

• This file / tutorial has been especially useful in figuring out how the OPLS file is organized: http://www.ks.uiuc.edu/Training/Tutorials/science/topology/topology-tutorial.pdf
• IC means internal coordinates
• RES = residue; PRES = patch residue (not sure what this means exactly. Maybe a residue that's added that's unknown in nature?)
    correction: It turns out that patch residues are "similar to regular residues, except that they are used to modify an existing residue. Any atoms they contain replace or add to those in the residue; they can also remove atoms." From here: https://integrativemodeling.org/2.3.0/doc/html/classIMP_1_1atom_1_1CHARMMPatch.html

The actual JCTC paper is very helpful in indicating what possible bond types there are.
    Paper: https://pubs.acs.org/doi/pdf/10.1021/acs.jctc.5b00356
    Supplement: https://pubs.acs.org/doi/suppl/10.1021/acs.jctc.5b00356/suppl_file/ct5b00356_si_001.pdf

There are 6 dihedral angles for OPLS-AA:
    phi:    C   -   N   -   CA  -   C
    psi:    N   -   CA  -   C   -   N
    phi':   C   -   N   -   CA  -   CB
    psi':   CB  -   CA  -   C   -   N
    phi'':  C   -   N   -   CA  -   HA
    psi'':  HA  -   CA  -   C   -   N


Something to keep in mind is that it may be good to eventually represent everything as a tree of some sort.