pandoc -o Casimir_honeycomb_drive.docx -t docx Casimir_honeycomb_drive.tex
# https://ja01.chem.buffalo.edu/etcetera/latex-pandoc-word.html

# pandoc -F pandoc-crossref -M autoEqnLabels -M tableEqns --citeproc -t markdown-citations  -s pandoc-template.tex -f latex -t docx -o converted_2.docx  --bibliography=thecitations.bib  --csl=journal-of-the-american-chemical-society.csl

# pandoc -F pandoc-crossref -M autoEqnLabels -M tableEqns --citeproc -t markdown-citations  -s Casimir_honeycomb_drive.tex -f latex -t docx -o converted.docx # --bibliography=thecitations.bib # --csl=journal-of-the-american-chemical-society.csl

# pandoc: unrecognized option `--citeproc'

# pandoc: pandoc-template.tex: openFile: does not exist (No such file or directory)

# bpandoc: Error running filter pandoc-crossref:
# Could not find executable 'pandoc-crossref'.

# https://github.com/lierdakil/pandoc-crossref
# https://github.com/lierdakil/pandoc-crossref/releases/tag/v0.3.17.0

# WARNING: pandoc-crossref was compiled with pandoc 3.1.8 but is being run through 1.19.2.4. This is not supported. Strange things may (and likely will) happen silently.
# pandoc-crossref: Error in $: Incompatible API versions: encoded with [1,17,0,5] but attempted to decode with [1,23,1].
# CallStack (from HasCallStack):
#   error, called at src/Text/Pandoc/JSON.hs:120:27 in pandoc-types-1.23.1-7RMcFb07Rsv1aOeScOW1Ry:Text.Pandoc.JSON
# pandoc: Error running filter pandoc-crossref
# Filter returned error status 1




# https://stackoverflow.com/questions/61100045/how-to-install-stable-and-fresh-pandoc-on-ubuntu

# Purge apt installed, you must "undo" the bad instruction.
# apt purge pandoc
# Get fresh .deb at Pandoc's git/releases. Example:
# wget https://github.com/jgm/pandoc/releases/download/2.9.2.1/pandoc-2.9.2.1-1-amd64.deb
# Install. Example:
# sudo dpkg -i pandoc-2.9.2.1-1-amd64.deb


# https://github.com/jgm/pandoc/releases
# https://github.com/jgm/pandoc/releases/tag/3.1.8

#wget https://github.com/jgm/pandoc/releases/download/3.1.8/pandoc-3.1.8-1-amd64.deb
#sudo dpkg -i pandoc-3.1.8-1-amd64.deb



#pandoc -F pandoc-crossref -M autoEqnLabels -M tableEqns -t markdown-citations -f latex -t docx -o converted.docx  --bibliography=thecitations.bib  --csl=journal-of-the-american-chemical-society.csl #Casimir_honeycomb_drive.tex
