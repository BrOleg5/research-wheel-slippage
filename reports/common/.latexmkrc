# source: https://stackoverflow.com/a/64860253
# add path to style and configuration files
ensure_path('TEXINPUTS', '../common//');

# source: https://github.com/nasa/nasa-latex-docs/blob/master/support/latexmk/latexmkrc
# set pdf viewer
if ($^O =~ /MSWin32/) {
    # Use Windows file associations
    $pdf_previewer = 'start %S';
} elsif ($^O =~ /linux/) {
    # Use default pdf viewer evince
    $pdf_previewer = 'start evince %O %S';
}

$xelatex = 'xelatex -synctex=1 -interaction=nonstopmode %O %S';
$pdf_mode = 5;
$dvi_mode = 0;
$postscript_mode = 0;
$preview_mode = 1;