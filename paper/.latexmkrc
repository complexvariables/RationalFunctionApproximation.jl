unshift @ARGV, '-pdf', '-bibtex';
sub build_header {
  system("ruby ./prep.rb")
}

build_header()
