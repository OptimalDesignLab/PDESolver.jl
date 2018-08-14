# Perl script for checking a PDESolver input dictionary for duplicates
# Only intended to print output, not do any fixing of the problem.
#
# Can be run outside of PDESolver with:
#   perl check_dict_for_dupes.pl dict.jl
# where dict.jl is some properly formatted PDESolver input dictionary

my $counter = 0;      # This counter should correspond to the number of input options 
                      #   specified in the input dictionary.
my %hash_for_check;   # We want to assign the dict key-value pairs from the file
                      #   to this hash. Perl hashes act like Julia dicts, so duplicates get
                      #   overwritten. We will compare these two lengths.

while (<>) {          # while there are lines of stdin (loop over file lines)
  no warnings;
  if (/Dict{Any,Any}/) {    # skip the first line
    next;
  }
  if (/^\s*\)\s*/) {        # skip the last line (only a parenthesis possibly surrounded by whitespace)
    next;
  }
  if (/^\s*$/) {            # skip blank lines (possibly containing whitespace)
    next;
  }
  if (/^\s*#/) {            # skip commented-out lines (line starts w/ possible whitespace followed by #)
    next;
  }
  # match:
  #   beginning of line
  #   possible whitespace
  #   quotation mark
  #   at least one character that isn't a quotation mark    -> stored in $1
  #   quotation mark
  #   possible whitespace
  #   the string "=>"
  #   possible whitespace
  #   at least one character, unspecified type    -> stored in $2
  #   possible whitespace
  #   comma
  #   possible whitespace
  #   end of line
  # The intent here is to match exactly the key & value pair of each line in the dictionary file.
  # This is purposefully fairly stringent.
  if (/^\s*\"([^"]+)\"\s*=>\s*(.+)\s*,\s*$/) {      
    # build a perl hash of all the key-value pairs (similar to julia Dict)
    $hash_for_check{$1} = $2;
  }
  $counter++;
}

# If the hash length is less than the number of key-value pairs,
#   then there probably is a duplicate.
#   This is because perl hashes and julia dicts are built in a similar fashion:
#   if a key already exists, it is just overwritten if it is assigned again.
if (length(%hash_for_check) < $counter) {
  print "=================================================\n";
  print "  Input dictionary: possible duplicate spotted.\n";
  print "=================================================\n";
}
