#!/usr/bin/perl

# ===========================================================================
#
#                            PUBLIC DOMAIN NOTICE
#            National Center for Biotechnology Information (NCBI)
#
#  This software/database is a "United States Government Work" under the
#  terms of the United States Copyright Act.  It was written as part of
#  the author's official duties as a United States Government employee and
#  thus cannot be copyrighted.  This software/database is freely available
#  to the public for use. The National Library of Medicine and the U.S.
#  Government do not place any restriction on its use or reproduction.
#  We would, however, appreciate having the NCBI and the author cited in
#  any work or product based on this material.
#
#  Although all reasonable efforts have been taken to ensure the accuracy
#  and reliability of the software and data, the NLM and the U.S.
#  Government do not and cannot warrant the performance or results that
#  may be obtained by using this software or data. The NLM and the U.S.
#  Government disclaim all warranties, express or implied, including
#  warranties of performance, merchantability or fitness for any particular
#  purpose.
#
# ===========================================================================
#
# File Name:  xtract
#
# Author:  Jonathan Kans
#
# Version Creation Date:   8/20/12
#
# ==========================================================================

# Entrez Direct - EDirect

# use strict;
use warnings;

# definitions

use constant false => 0;
use constant true  => 1;

# xtract version number

$version = "2.9999";

# initialize memory hash

my %memory = ();

# initialize synopsis variables

my %synopsis_acc = ();
my $synopsis_max = 0;

# initial level is global

my $initial_level = 7;

# utility subroutines

sub convert_http {
  my $str = shift (@_);
  $str =~ s/\&amp\;/\&/g;
  $str =~ s/\&apos\;/\'/g;
  $str =~ s/\&gt\;/\>/g;
  $str =~ s/\&lt\;/\</g;
  $str =~ s/\&quot\;/\"/g;
  return $str;
}

sub convert_slash {
  my $str = shift (@_);
  $str =~ s/\\t/\t/g;
  $str =~ s/\\n/\n/g;
  $str =~ s/\\r/\n/g;
  return $str;
}

# xtract parses and extracts data from XML, driven by command-line arguments

sub process_element {

  my $line = shift (@_);
  my $val = shift (@_);
  my $prev = shift (@_);
  my $pfx = shift (@_);
  my $sfx = shift (@_);
  my $sep = shift (@_);
  my $cmmd = shift (@_);
  my $variable = shift (@_);
  my $indx = shift (@_);
  my $sclr = shift (@_);

  if ( $val eq "" ) {
    return false;
  }

  # -element "*" writes current XML object

  if ( $val eq "*" ) {
    print "$prev";
    print "$pfx";
    print "$line";
    print "$sfx";
    return true;
  }

  my @accum = ();

  # commas are used to select unrelated items to show between pfx and sfx

  my @itms = split (',', $val);
  foreach $itm (@itms) {

    my @atts = ();
    my @vals = ();
    my @working = ();

    my $do_count = false;
    my $do_length = false;


    # look for # or % modifiers before XML element specification

    if ( $itm =~ /^#(.+)/ ) {

      # items in quotations starting with # sign display count of matches

      $do_count = true;
      $itm = $1;

    } elsif ( $itm =~ /^%(.+)/ ) {

      # items in quotations starting with % sign display sum of matched string lengths

      $do_length = true;
      $itm = $1;
    }


    # INSDSeq special cases - &Feature and &Qualifier

    if ( $itm =~ /^\&Feature$/ ) {

      $itm = "INSDFeature_key";

    } elsif ( $itm =~ /^\&Qualifier$/ ) {

      $itm = "INSDQualifier_value";
    }


    # -element pattern recognition

    if ( $itm =~ /^&([A-Z0-9]+)$/ ) {

      # match stored &VARIABLE

      my $idx = $1;
      if ( exists ($memory{$idx}) and $memory{"$idx"} ne "" ) {
        my $val = $memory{$idx};
        $dec = convert_http($val);
        push (@working, $dec);
      }

    } elsif ( $itm =~ /^\?([A-Z]+)\((.*)\)/ ) {

      # match exploration flag (alias to print is in parentheses) (undocumented)

      my $itm = $1;
      my $rep = $2;

      $dec = convert_http($rep);

      # first and last elements are always present, and are the same for a one-element list

      if ( $itm eq "FIRST" and $indx == 1 ) {
        push (@working, $dec);
      } elsif ( $itm eq "LAST" and $indx == $sclr ) {
        push (@working, $dec);

      # single, head, inner, tail are mutually exclusive situations

      } elsif ( $itm eq "SINGLE" and $sclr == 1 ) {
        push (@working, $dec);

      } elsif ( $itm eq "HEAD" and $indx == 1 and $sclr > 1 ) {
        push (@working, $dec);
      } elsif ( $itm eq "INNER" and $indx > 1 and $indx < $sclr ) {
        push (@working, $dec);
      } elsif ( $itm eq "TAIL" and $indx == $sclr and $sclr > 1 ) {
        push (@working, $dec);

      # even and odd are based only on index value, not list length

      } elsif ( $itm eq "EVEN" and ( $indx % 2) == 0 ) {
        push (@working, $dec);
      } elsif ( $itm eq "ODD" and ( $indx % 2) == 1 ) {
        push (@working, $dec);
      }

    } elsif ( $itm =~ /^@(.+)/ ) {

      # match unqualified @attribute

      my $att = $1;
      @atts = ($line =~ /<\S+ ([^>]+)>/g);
      foreach $val (@atts) {
        if ( $val =~ /$att=\"([^\"]+)\"/ ) {
          $dec = convert_http($1);
          push (@working, $dec);
        }
      }

    } elsif ( $itm =~ /(.+)@(.+)/ ) {

      # match qualified tag@attribute

      my $tag = $1;
      my $att = $2;
      @atts = ($line =~ /<$tag ([^>]+)>/g);
      foreach $val (@atts) {
        if ( $val =~ /$att=\"([^\"]+)\"/ ) {
          $dec = convert_http($1);
          push (@working, $dec);
        }
      }

    } elsif ( $itm =~ /(.+)\/(.+)/ ) {

      # match parent/child by tracking depth level of tokens

      my $tag = $1;
      my $chd = $2;
      @vals = ($line =~ /<$tag(?:\s+.+?)?>(.+?)<\/$tag>/g);
      foreach $val (@vals) {
        my $lvl = 0;
        my @tokens = split (/(?<=>)(?=<)/, $val);
        foreach $tkn (@tokens) {
          if ( $lvl < 0 ) {
            # ignore remainder if level drops below zero
          } elsif ( $tkn =~ /^<.+?>.+?<\/.+?>$/ ) {
            # content-containing tag does not change level
            if ( $tkn =~ /<$chd(?:\s+.+?)?>(.+?)<\/$chd>/ ) {
              # child tag matches, only accept if immediately below parent
              if ( $lvl == 0 ) {
                $dec = convert_http($1);
                push (@working, $dec);
              }
            }
          } elsif ( $tkn =~ /\/>$/ ) {
            # self-closing token does not change level
          } elsif ( $tkn =~ /^<[^\/]/ ) {
            # open tag increments level
            $lvl++;
          } elsif ( $tkn =~ /^<\// ) {
            # close tag decrements level
            $lvl--;
          }
        }
      }

    } elsif ( $itm =~ /(.+)\((.*)\)/ ) {

      # match self-closing tag (alias to print is in parentheses)

      my $tag = $1;
      my $rep = $2;
      if ( $rep eq "" ) {
        $rep = $tag;
      }

      # [ is it possible to look for both forms in one pass ? ]

      @vals = ($line =~ /(<$tag.*?><\/$tag>)/g);

      if ( scalar @vals != 0 ) {

        # found alternative form ( <tag></tag> ) of self-closing tag

        foreach $val (@vals) {
          $dec = convert_http($rep);
          push (@working, $dec);
        }

      } else {

        # look for compact form ( <tag/> ) if no match to alternative form

        if ( scalar @vals == 0 ) {
          @vals = ($line =~ /(<$tag.*?\/>)/g);
          foreach $val (@vals) {
            $dec = convert_http($rep);
            push (@working, $dec);
          }
        }
      }

    } else {

      # match tag with contents

      my $tag = $itm;
      @vals = ($line =~ /<$tag(?:\s+.+?)?>(.+?)<\/$tag>/g);
      foreach $val (@vals) {
        $dec = convert_http($val);
        push (@working, $dec);
      }
    }


    my $num = scalar @working;

    # if -first, -last, etc., remove all but matching subset

    if ( $num < 1 ) {

      # no matches collected, skip

    } elsif ( $cmmd eq "-first" ) {

      @working = splice (@working, 0, 1);
      $num = scalar @working;

    } elsif ( $cmmd eq "-last" ) {

      @working = splice (@working, $num - 1, 1);
      $num = scalar @working;

    } elsif ( $cmmd eq "-single" ) {

      if ( $num > 1 ) {
        $num = 0;
      }

    } elsif ( $cmmd eq "-head" ) {

      if ( $num < 2 ) {
        $num = 0;
      } else {
        @working = splice (@working, 0, 1);
        $num = scalar @working;
      }

    } elsif ( $cmmd eq "-inner" ) {

      if ( $num < 3 ) {
        $num = 0;
      } else {
        @working = splice (@working, 1, $num - 2);
        $num = scalar @working;
      }

    } elsif ( $cmmd eq "-tail" ) {

      if ( $num < 2 ) {
        $num = 0;
      } else {
        @working = splice (@working, $num - 1, 1);
        $num = scalar @working;
      }
    }


    # store current results into accumulator

    if ( $num == 0 ) {

      # no printable results

    } elsif ( $do_count ) {

      push (@accum, $num);

    } elsif ( $do_length ) {

      my $len = 0;
      foreach $dec (@working) {
        $len += length ( $dec );
      }

      push (@accum, $len);

    } elsif ( $cmmd eq "-even" ) {

      my $take = false;
      foreach $dec (@working) {
        if ( $take ) {
          push (@accum, $dec);
        }
        $take = (! $take);
      }

    } elsif ( $cmmd eq "-odd" ) {

      my $take = true;
      foreach $dec (@working) {
        if ( $take ) {
          push (@accum, $dec);
        }
        $take = (! $take);
      }

    } else {

      foreach $dec (@working) {
        push (@accum, $dec);
      }
    }
  }


  # process results

  if ( scalar @accum == 0 ) {

    return false;

  } elsif ( $variable ne "") {

    # store accumulated components into variable

    my $str = $pfx;
    my $between = "";
    foreach $dec (@accum) {
      $str .= "$between";
      $between = $sep;
      $str .= "$dec";
    }
    $str .= "$sfx";

    $memory{"$variable"} = "$str";

  } else {

    # format accumulated elements

    my $before = $prev . $pfx;
    my $between = "";
    my $after = "";
    foreach $dec (@accum) {
      print "$before";
      $before = "";
      print "$between";
      $between = $sep;
      $after = $sfx;
      print "$dec";
    }
    print "$after";
    return true;
  }

  return false;
}

sub process_flags {

  my $line = shift (@_);
  my @args = @{shift (@_)};
  my $tab = shift (@_);
  my $ret = shift (@_);
  my $indx = shift (@_);
  my $sclr = shift (@_);

  my $pfx = "";
  my $sfx = "";
  my $sep = "\t";

  my $col = "\t";
  my $lin = "\n";

  my $okay = false;
  my $cmmd = "";

  my $max = scalar @args;

  for ( my $i = 0; $i < $max; $i++ ) {

    # read command-line arguments

    my $val = $args[$i];

    # if new argument starts with hyphen

    if ( $val =~ /^\-/ ) {

      # not an -element value, so set $okay flag to false

      $okay = false;

      # all flags should be followed by at least one value

      if ( $i + 1 >= $max and $val ne "-outline" and $val ne "-synopsis" ) {

        # warn that the last flag has no value

        print STDERR "No argument after '$val'\n";

        # break out of for loop to return from function

        last;

        # otherwise safe to increment with $i++ and then dereference $args [$i]
      }
    }

    # print outline of XML structure

    if ( $val eq "-outline" ) {

      my $lvl = 0;
      my $dpth = 0;

      my @tokens = split (/(?<=>)(?=<)/, $line);
      foreach $tkn (@tokens) {
        if ( $lvl < 0 ) {
          # ignore remainder if level drops below zero
        } elsif ( $tkn =~ /^<.+?>.+?<\/(.+?)>$/ ) {
          # content-containing tag
          if ( $tkn =~ /^<(\S+).*?>.*>$/ ) {
            $tkn = $1;
            $lvl++;
            $dpth++;
            for ( $i = 0; $i < $dpth; $i++ ) {
              print "  ";
            }
            print "$tkn\n";
            $lvl--;
            $dpth--;
          }
        } elsif ( $tkn =~ /^<.+?\/>$/ ) {
          # self-closing token
          if ( $tkn =~ /^<(\S+).*?\/>$/ ) {
            $tkn = $1;
            $lvl++;
            $dpth++;
            for ( $i = 0; $i < $dpth; $i++ ) {
              print "  ";
            }
            print "$tkn\n";
            $lvl--;
            $dpth--;
          }
        } elsif ( $tkn =~ /^<[^\/]/ ) {
          # open tag increments level
          $lvl++;
          $dpth++;
          if ( $tkn =~ /^<(\S+).*?>$/ ) {
            $tkn = $1;
            if ( $tkn ne "?xml" and
                 $tkn ne "!DOCTYPE" and
                 $tkn ne "eSummaryResult" and
                 $tkn ne "eLinkResult" and
                 $tkn ne "eInfoResult" and
                 $tkn ne "PubmedArticleSet" and
                 $tkn ne "DocumentSummarySet" and
                 $tkn ne "INSDSet" and
                 $tkn ne "Entrezgene-Set" and
                 $tkn ne "TaxaSet" ) {
              for ( $i = 0; $i < $dpth; $i++ ) {
                print "  ";
              }
              print "$tkn\n";
            } else {
              $dpth--;
            }
          }
        } elsif ( $tkn =~ /^<\// ) {
          # close tag decrements level
          $lvl--;
          $dpth--;
        }
      }

    # collect summary of XML paths, print at end

    } elsif ( $val eq "-synopsis" ) {

      my @arr = ();
      my $lvl = 0;

      my @tokens = split (/(?<=>)(?=<)/, $line);
      foreach $tkn (@tokens) {
        my $record_level = 0;
        if ( $lvl < 0 ) {
          # ignore remainder if level drops below zero
        } elsif ( $tkn =~ /^<.+?>.+?<\/(.+?)>$/ ) {
          # content-containing tag
          if ( $tkn =~ /^<(\S+).*?>.*>$/ ) {
            $tkn = $1;
            $lvl++;
            $arr[$lvl] = $tkn;
            $record_level = $lvl;
            $lvl--;
          }
        } elsif ( $tkn =~ /^<.+?\/>$/ ) {
          # self-closing token
          if ( $tkn =~ /^<(\S+).*?\/>$/ ) {
            $tkn = $1;
            $lvl++;
            $arr[$lvl] = $tkn;
            $record_level = $lvl;
            $lvl--;
          }
        } elsif ( $tkn =~ /^<[^\/]/ ) {
          # open tag increments level
          $lvl++;
          if ( $tkn =~ /^<(\S+).*?>$/ ) {
            $tkn = $1;
            $arr[$lvl] = $tkn;
            $record_level = $lvl;
          }
        } elsif ( $tkn =~ /^<\// ) {
          # close tag decrements level
          $lvl--;
        }
        if ( $record_level > 0 ) {
          my $str = "";
          my $pfx = "";
          for ( $i = 1; $i <= $record_level; $i++ ) {
            $itm = $arr[$i];
            if ( $itm ne "?xml" and
                 $itm ne "!DOCTYPE" and
                 $itm ne "eSummaryResult" and
                 $itm ne "eLinkResult" and
                 $itm ne "eInfoResult" and
                 $itm ne "PubmedArticleSet" and
                 $itm ne "DocumentSummarySet" and
                 $itm ne "INSDSet" and
                 $itm ne "Entrezgene-Set" and
                 $itm ne "TaxaSet" ) {
              $str .= "$pfx";
              $pfx = "/";
              $str .= "$arr[$i]";
            }
          }
          if ( $str ne "" ) {
            my $val = 0;
            if ( exists ($synopsis_acc{"$str"}) ) {
              $val = $synopsis_acc{"$str"};
            }
            $val++;
            $synopsis_acc{"$str"} = $val;
            if ( $val > $synopsis_max ) {
              $synopsis_max = $val;
            }
            print "$str\n"
          }
        }
      }

    # look for arguments that adjust output formatting

    } elsif ( $val eq "-pfx" ) {
      $i++;
      $pfx = convert_slash ( $args[$i] );
    } elsif ( $val eq "-sfx" ) {
      $i++;
      $sfx = convert_slash ( $args[$i] );
    } elsif ( $val eq "-sep" ) {
      $i++;
      $sep = convert_slash ( $args[$i] );
    } elsif ( $val eq "-tab" ) {
      $i++;
      $col = convert_slash ( $args[$i] );
    } elsif ( $val eq "-ret" ) {
      $i++;
      $lin = convert_slash ( $args[$i] );
    } elsif ( $val eq "-lbl" ) {
      $i++;
      print "$tab";
      $lbl = convert_slash ( $args[$i] );
      print "$lbl";
      $tab = $col;
      $ret = $lin;

    # look for -element or limiting variants, non-hyphenated arguments to follow

    } elsif ( $val eq "-element" or
              $val eq "-first" or
              $val eq "-last" or
              $val eq "-even" or
              $val eq "-odd" or
              $val eq "-single" or
              $val eq "-head" or
              $val eq "-inner" or
              $val eq "-tail" ) {

      $okay = true;
      $cmmd = $val;

    # hyphen followed by capital letters or digits is a variable

    } elsif ( $val =~ /^\-([A-Z0-9]+)([a-z]?)$/ ) {

      my $variable = $1;
      my $first_or_last = $2;

      $cmmd = "-element";

      # detect f, l, etc., variable suffix for -first, -last, etc. (undocumented)

      if ( $first_or_last eq "f" ) {
        $cmmd = "-first";
      } elsif ( $first_or_last eq "l" ) {
        $cmmd = "-last";
      } elsif ( $first_or_last eq "e" ) {
        $cmmd = "-even";
      } elsif ( $first_or_last eq "o" ) {
        $cmmd = "-odd";
      } elsif ( $first_or_last eq "s" ) {
        $cmmd = "-single";
      } elsif ( $first_or_last eq "h" ) {
        $cmmd = "-head";
      } elsif ( $first_or_last eq "i" ) {
        $cmmd = "-inner";
      } elsif ( $first_or_last eq "t" ) {
        $cmmd = "-tail";
      }

      # clear the variable to avoid old value persisting if no subsequent match
      # commented out in favor of clearing by assigning a literal value containing an empty string

      # $memory{"$variable"} = "";

      $i++;
      $val = $args[$i];

      # set variable to new value

      if ( $val =~ /^\((.*)\)$/ ) {

        # parentheses contain literal value for variable

        $val = convert_slash ( $1 );
        $memory{"$variable"} = "$val";

      } else {

        # if XML tag name, record what would otherwise be printed by -element

        process_element ( $line, $val, $tab, $pfx, $sfx, $sep, $cmmd, $variable, $indx, $sclr );

      }

    # reality check on unsupported arguments

    } elsif ( $val =~ /^\-/ ) {

      print STDERR "Unrecognized argument '$val'\n";

    } elsif ( $okay ) {

      # process -element values

      if ( process_element ( $line, $val, $tab, $pfx, $sfx, $sep, $cmmd, "", $indx, $sclr )) {
        $tab = $col;
        $ret = $lin;
      }

    } else {

      print STDERR "No -element before '$val'\n";

    }
  }

  return $tab, $ret;
}

sub has_content {

  my $line = shift (@_);
  my $str = shift (@_);
  my $indx = shift (@_);
  my $sclr = shift (@_);

  if ( $str eq "" or $str eq "&" ) {
    return false;
  }


  # match by non-empty variable value (undocumented)

  if ( $str =~ /^&([A-Z0-9]+)$/ ) {
    $str = $1;
    if ( exists ($memory{"$str"}) and $memory{"$str"} ne "" ) {
      return true;
    }
    return false;
  }

  # remove leading backslash used to prevent content from being mistaken for variable

  $str =~ s/^\\&/&/;


  # match by exploration flag (undocumented)

  if ( $str =~ /^\?([A-Z]+)/ ) {
    $str = $1;

    # first and last elements are always present, and are the same for a one-element list

    if ( $str eq "FIRST" and $indx == 1 ) {
      return true;
    } elsif ( $str eq "LAST" and $indx == $sclr ) {
      return true;

    # single, head, inner, tail are mutually exclusive situations

    } elsif ( $str eq "SINGLE" and $sclr == 1 ) {
      return true;

    } elsif ( $str eq "HEAD" and $indx == 1 and $sclr > 1 ) {
      return true;
    } elsif ( $str eq "INNER" and $indx > 1 and $indx < $sclr ) {
      return true;
    } elsif ( $str eq "TAIL" and $indx == $sclr and $sclr > 1 ) {
      return true;

    # even and odd are based only on index value, not list length

    } elsif ( $str eq "EVEN" and ( $indx % 2) == 0 ) {
      return true;
    } elsif ( $str eq "ODD" and ( $indx % 2) == 1 ) {
      return true;
    }
  }

  # remove leading backslash used to prevent content from being mistaken for flag

  $str =~ s/^\\\?/?/;


  # exact match by data contents, not using regular expression, case insensitive

  if ( index ( uc($line), uc($str) ) >= 0 ) {
    return true;
  }

  return false;
}

sub has_element {

  my $line = shift (@_);
  my $str = shift (@_);

  if ( $str eq "" ) {
    return false;
  }


  # track depth level of tokens if pattern was parent/child (undocumented)
  
  if ( $str =~ /(.+)\/(.+)/ ) {
    my $ptrn = $1;
    my $chld = $2;

    my @vals = ();
    my $lvl = 0;
    my $dpth = 0;
    my $in_pat = false;

    my @tokens = split (/(?<=>)(?=<)/, $line);
    foreach $tkn (@tokens) {

      if ( $lvl < 0 ) {
        # ignore remainder if level drops below zero
      } elsif ( $tkn =~ /^<.+?>.+?<\/.+?>$/ ) {
        # content-containing tag does not change level
        if ( $in_pat and $lvl == $dpth ) {
          push (@vals, $tkn);
        }
      } elsif ( $tkn =~ /\/>$/ ) {
        # self-closing token does not change level
        if ( $in_pat and $lvl == $dpth ) {
          push (@vals, $tkn);
        }
      } elsif ( $tkn =~ /^<[^\/]/ ) {
        # open tag increments level
        $lvl++;
        if ( $tkn =~ /<$ptrn(?:\s+.+?)?>/ and $dpth == 0 ) {
          $in_pat = true;
          $dpth = $lvl;
          push (@vals, $tkn);
        }
      } elsif ( $tkn =~ /^<\// ) {
        # close tag decrements level
        if ( $in_pat and $tkn =~ /<\/$ptrn>/ and $lvl == $dpth ) {
          push (@vals, $tkn);
          $in_pat = false;
          $dpth = 0;
        }
        $lvl--;
      }
    }

    if ( scalar @vals > 0 ) {
      my $sub = join ("", @vals);
      if ( has_element ( $sub, $chld )) {
        return true;
      }
    }

    return false;
  }


  # match by specific variable value (undocumented)

  if ( $str =~ /^\&([A-Z0-9]+)\:(.+)$/ ) {
    my $ptrn = $1;
    my $chld = $2;

    if ( exists ($memory{"$ptrn"}) and $memory{"$ptrn"} eq "$chld" ) {
      return true;
    }

    return false;
  }


  # match by non-empty variable value (undocumented)

  if ( $str =~ /^\&([A-Z0-9]+)$/ ) {
    $str = $1;

    if ( exists ($memory{"$str"}) and $memory{"$str"} ne "" ) {
      return true;
    }

    return false;
  }


  # regular expression search for element@attribute:value

  if ( $str =~ /(.+)@(.+)\:(.+)/ ) {
    my $ptrn = $1;
    my $attr = $2;
    my $chld = $3;

    if ( $line =~ /<$ptrn ([^>]+)>/i ) {
      my $atts = $1;
      if ( $atts =~ /$attr=\"$chld\"/ ) {
        return true;
      }
    }

    return false;
  }


  # regular expression search for @attribute:value

  if ( $str =~ /@(.+)\:(.+)/ ) {
    my $attr = $1;
    my $chld = $2;

    if ( $line =~ /<\S+ ([^>]+)>/i ) {
      my $atts = $1;
      if ( $atts =~ /$attr=\"$chld\"/ ) {
        return true;
      }
    }

    return false;
  }


  # regular expression search for element@attribute

  if ( $str =~ /(.+)@(.+)/ ) {
    my $ptrn = $1;
    my $attr = $2;

    if ( $line =~ /<$ptrn ([^>]+)>/i ) {
      my $atts = $1;
      if ( $atts =~ /$attr=\"([^\"]+)\"/ ) {
        return true;
      }
    }

    return false;
  }


  # regular expression search for @attribute

  if ( $str =~ /@(.+)/ ) {
    my $attr = $1;

    if ( $line =~ /<\S+ ([^>]+)>/i ) {
      my $atts = $1;
      if ( $atts =~ /$attr=\"([^\"]+)\"/ ) {
        return true;
      }
    }

    return false;
  }


  # INSDSeq special cases - &Feature:name and &Qualifier:name

  if ( $str =~ /^\&Feature\:(.+)/ ) {

    $str = "INSDFeature_key:" . "$1";

  } elsif ( $str =~ /^\&Qualifier\:(.+)/ ) {

    $str = "INSDQualifier_name:" . "$1";
  }


  # regular expression search for element:value

  if ( $str =~ /(.+)\:(.+)/ ) {
    my $ptrn = $1;
    my $chld = $2;

    if ( $line =~ /<$ptrn(?:\s+.+?)?>$chld<\/$ptrn>/i ) {
      return true;
    }

    return false;
  }


  # regular expression search for element tag

  if ( $line =~ /<$str(?:\s+.+?)?>/i ) {
    return true;
  }

  # also test self-closing tag

  if ( $line =~ /<$str(?:\s+.+?)?\/>/i ) {
    return true;
  }

  return false;
}

sub process_level {

  my $line = shift (@_);
  my $par = shift (@_);
  my $chd = shift (@_);
  my @args = @{shift (@_)};
  my $tab = shift (@_);
  my $ret = shift (@_);
  my $level = shift (@_);
  my $indx = shift (@_);
  my $sclr = shift (@_);


  # handle heterogeneous child blocks if pattern was parent/*

  if ( $chd eq "*" and $par ne "" and $par ne "*" ) {

    my @vals = ();
    my $lvl = 0;
    my $in_pat = false;
    my $pfx = "";
    my @working = ();

    if ( $line =~ /^<$par(?:\s+.+?)?>(.+?)<\/$par>$/ ) {
      $line = $1;
    }

    my @tokens = split (/(?<=>)(?=<)/, $line);
    foreach $tkn (@tokens) {

      if ( $lvl < 0 ) {
        # ignore remainder if level drops below zero
      } elsif ( $tkn =~ /^<.+?>.+?<\/.+?>$/ ) {
        # content-containing tag does not change level
        if ( $lvl > 0 ) {
          push (@working, $tkn);
        } else {
          # match to full object, push contents directly to final array
          push (@vals, $tkn);
        }
      } elsif ( $tkn =~ /\/>$/ ) {
        # self-closing token does not change level
        if ( $lvl > 0 ) {
          push (@working, $tkn);
        } else {
          push (@vals, $tkn);
        }
      } elsif ( $tkn =~ /^<[^\/]/ ) {
        # open tag increments level
        push (@working, $tkn);
        $lvl++;
      } elsif ( $tkn =~ /^<\// ) {
        # close tag decrements level
        push (@working, $tkn);
        $lvl--;
        if ( $lvl == 0 ) {
          my $str = join ("", @working);
          push (@vals, $str);
          @working = ();
        }
      }
    }
    if ( scalar @working > 0 ) {
      my $str = join ("", @working);
      push (@vals, $str);
      @working = ();
    }

    foreach $val (@vals) {
      if ( $level == $initial_level ) {
        print "$pfx";
        %memory = ();
      }
      ( $tab, $ret ) = process_level ( $val, "", "", \@args, $tab, $ret, $level, $indx, $sclr );
      if ( $level == $initial_level ) {
        $pfx = $ret;
        $tab = "";
      }
    }

    return $tab, $ret;
  }

  # track depth level of tokens if pattern was parent/child

  if ( $chd ne "" ) {

    my @vals = ();
    my $lvl = 0;
    my $dpth = 0;
    my $in_pat = false;
    my $pfx = "";
    my @working = ();

    my @tokens = split (/(?<=>)(?=<)/, $line);
    foreach $tkn (@tokens) {

      if ( $lvl < 0 ) {
        # ignore remainder if level drops below zero
      } elsif ( $tkn =~ /^<.+?>.+?<\/.+?>$/ ) {
        # content-containing tag does not change level
        if ( $in_pat ) {
          push (@working, $tkn);
        } elsif ( $tkn =~ /<$chd(?:\s+.+?)?>(.+?)<\/$chd>/ and $lvl < 2 ) {
          # match to full object, push directly to final array
          push (@vals, $tkn);
        }
      } elsif ( $tkn =~ /\/>$/ ) {
        # self-closing token does not change level
        if ( $in_pat ) {
          push (@working, $tkn);
        }
      } elsif ( $tkn =~ /^<[^\/]/ ) {
        # open tag increments level
        if ( $in_pat ) {
          push (@working, $tkn);
        } elsif ( $tkn =~ /<$chd(?:\s+.+?)?>/ and $lvl < 2 ) {
          $in_pat = true;
          $dpth = $lvl;
          push (@working, $tkn);
        }
        $lvl++;
      } elsif ( $tkn =~ /^<\// ) {
        # close tag decrements level
        $lvl--;
        if ( $in_pat ) {
          push (@working, $tkn);
        }
        if ( $tkn =~ /<\/$chd>/ ) {
          if ( $lvl <= $dpth ) {
            $in_pat = false;
            my $str = join ("", @working);
            push (@vals, $str);
            @working = ();
          }
        }
      }
    }
    if ( scalar @working > 0 ) {
      my $str = join ("", @working);
      push (@vals, $str);
      @working = ();
    }

    foreach $val (@vals) {
      if ( $level == $initial_level ) {
        print "$pfx";
        %memory = ();
      }
      ( $tab, $ret ) = process_level ( $val, "", "", \@args, $tab, $ret, $level, $indx, $sclr );
      if ( $level == $initial_level ) {
        $pfx = $ret;
        $tab = "";
      }
    }

    return $tab, $ret;
  }


  my $max = scalar @args;


  # allow conditional arguments in any order

  my $go_on = true;
  while ( $go_on ) {

    if ( $max > 2 and $args[0] eq "-position" ) {

      # -position command filters by object position in list

      my $pstn = $args[1];
      $max -= 2;
      @args = splice (@args, 2, $max);

      # first and last elements are always present, and are the same for a one-element list

      if ( $pstn eq "first" and  $indx > 1 ) {
        return $tab, $ret;
      } elsif ( $pstn eq "last" and $indx < $sclr ) {
        return $tab, $ret;

      # single, head, inner, tail are mutually exclusive situations

      } elsif ( $pstn eq "single" and $sclr > 1 ) {
        return $tab, $ret;

      } elsif ( $pstn eq "head" and ( $indx > 1 or $sclr == 1 ) ) {
        return $tab, $ret;
      } elsif ( $pstn eq "inner" and ( $indx == 1 or $indx == $sclr ) ) {
        return $tab, $ret;
      } elsif ( $pstn eq "tail" and ( $indx < $sclr or $sclr == 1 ) ) {
        return $tab, $ret;

      # single and multiple are mutually exclusive situations

      } elsif ( $pstn eq "multiple" and $sclr == 1 ) {
        return $tab, $ret;

      # even and odd are based only on index value, not list length

      } elsif ( $pstn eq "even" and ( $indx % 2) == 1 ) {
        return $tab, $ret;
      } elsif ( $pstn eq "odd" and ( $indx % 2) == 0 ) {
        return $tab, $ret;

      # value is actual numeric index

      } elsif ( $pstn =~ /^[0-9]+$/ and $indx != $pstn ) {
        return $tab, $ret;
      }

    } elsif ( $max > 2 and $args[0] eq "-match" ) {

      # -match ... -and/-or ... command filters for indicated element tag [:value]

      my $required = 1;
      my $observed = 0;
      my $prevbool = "";

      my $mtch = $args[1];
      $max -= 2;
      @args = splice (@args, 2, $max);

      if ( has_element ( $line, $mtch )) {
        $observed++;
      }

      # "or" succeeds on any match, "and" requires all tests to match

      while ( $max > 2 and ( $args[0] eq "-or" or $args[0] eq "-and" )) {

        if ( $prevbool ne "" and $prevbool ne $args[0] ) {
          print STDERR "A mixture of -and and -or commands cannot follow -match\n";
          return $tab, $ret;
        }

        $prevbool = $args[0];
        if ( $prevbool eq "-and" ) {
          $required++;
        }

        $mtch = $args[1];
        $max -= 2;
        @args = splice (@args, 2, $max);

        if ( has_element ( $line, $mtch )) {
          $observed++;
        }
      }

      if ( $observed < $required ) {
        return $tab, $ret;
      }

    } elsif ( $max > 2 and $args[0] eq "-avoid" ) {

      # -avoid ... -and ... command filters against indicated element tag [:value]

      my $skip = $args[1];
      $max -= 2;
      @args = splice (@args, 2, $max);

      if ( has_element ( $line, $skip )) {
        return $tab, $ret;
      }

      # colloquial "and" is really logical "or" - any match bails out of function

      while ( $max > 2 and $args[0] eq "-and" ) {

        $skip = $args[1];
        $max -= 2;
        @args = splice (@args, 2, $max);

        if ( has_element ( $line, $skip )) {
          return $tab, $ret;
        }
      }

      if ( $max > 1 and $args[0] eq "-or" ) {
        print STDERR "The -or command cannot follow -avoid\n";
        return $tab, $ret;
      }

    } elsif ( $max > 2 and $args[0] eq "-present" ) {

      # -present ... -and/-or ... command filters for indicated data content

      my $required = 1;
      my $observed = 0;
      my $prevbool = "";

      my $mtch = $args[1];
      $max -= 2;
      @args = splice (@args, 2, $max);

      if ( has_content ( $line, $mtch, $indx, $sclr )) {
        $observed++;
      }

      # "or" succeeds on any match, "and" requires all tests to match

      while ( $max > 2 and ( $args[0] eq "-or" or $args[0] eq "-and" )) {

        if ( $prevbool ne "" and $prevbool ne $args[0] ) {
          print STDERR "A mixture of -and and -or commands cannot follow -present\n";
          return $tab, $ret;
        }

        $prevbool = $args[0];
        if ( $prevbool eq "-and" ) {
          $required++;
        }

        $mtch = $args[1];
        $max -= 2;
        @args = splice (@args, 2, $max);

        if ( has_content ( $line, $mtch, $indx, $sclr )) {
          $observed++;
        }
      }

      if ( $observed < $required ) {
        return $tab, $ret;
      }

    } elsif ( $max > 2 and $args[0] eq "-absent" ) {

      # -absent ... -and ... command filters against indicated data content

      my $skip = $args[1];
      $max -= 2;
      @args = splice (@args, 2, $max);

      if ( has_content ( $line, $skip, $indx, $sclr )) {
        return $tab, $ret;
      }

      # colloquial "and" is really logical "or" - any match bails out of function

      while ( $max > 2 and $args[0] eq "-and" ) {

        $skip = $args[1];
        $max -= 2;
        @args = splice (@args, 2, $max);

        if ( has_content ( $line, $skip, $indx, $sclr )) {
          return $tab, $ret;
        }
      }

      if ( $max > 1 and $args[0] eq "-or" ) {
        print STDERR "The -or command cannot follow -absent\n";
        return $tab, $ret;
      }

    } elsif ( $max > 1 and $args[0] eq "-trim" ) {

      # -trim removes first XML tag

      while ( $max > 1 and $args[0] eq "-trim" ) {

        $max--;
        @args = splice (@args, 1, $max);

        if ( $line =~ /^<.+?>(.+)$/ ) {
          $line = $1;
        }
      }


    } else {

    # no other conditional tests or trim commands, break loop

      $go_on = false;
    }
  }


  # at level 0 break recursion for nested organizers, process -element

  if ( $level < 1 ) {

    ( $tab, $ret ) = process_flags ( $line, \@args, $tab, $ret, $indx, $sclr );
    return $tab, $ret;
  }


  # get name of command for current organizer level

  # initial capital signifies tokenized search matching nesting levels

  # (tokenizing is significantly slower than regular expression pattern matching)

  my $name = "";
  my $capname = "";

  if ( $level == 7 ) {
    $name = "-division";
    $capname = "-Division";
  } elsif ( $level == 6) {
    $name = "-group";
    $capname = "-Group";
  } elsif ( $level == 5) {
    $name = "-branch";
    $capname = "-Branch";
  } elsif ( $level == 4) {
    $name = "-block";
    $capname = "-Block";
  } elsif ( $level == 3) {
    $name = "-section";
    $capname = "-Section";
  } elsif ( $level == 2) {
    $name = "-subset";
    $capname = "-Subset";
  } elsif ( $level == 1) {
    $name = "-unit";
    $capname = "-Unit";
  }


  # find chains of arguments between indicated command names

  my $start = 0;
  my $stop = 0;

  for ( $start = 0; $start < $max; $start = $stop ) {
    $stop = $start + 1;
    while ( $stop < $max and $args[$stop] ne $name and $args[$stop] ne $capname ) {
      $stop++;
    }

    # collect arguments between -group, -block, or -subset directives

    my @tmp = @args;

    if ( $tmp[$start] eq $name and $stop - $start >= 2 ) {

      my $pat = $tmp[$start + 1];
      my $chd = "";

      if ( $pat =~ /\*\/(.+)/ ) {
        # -block */Taxon replaces -trim
        $pat = $1;
        if ( $line =~ /^<.+?>(.+)$/ ) {
          $line = $1;
        }
      }

      # INSDSeq special cases - &Feature and &Qualifier

      if ( $pat =~ /^\&Feature$/ ) {
        $pat = "INSDFeature";
      } elsif ( $pat =~ /^\&Qualifier$/ ) {
        $pat = "INSDQualifier";
      }

      if ( $pat =~ /(.+)\/(.+)/ ) {
        $pat = $1;
        $chd = $2;
      }
      $start += 2;
      my @group = splice (@tmp, $start, $stop - $start);

      # normal non-greedy pattern matching for non-nested XML

      my @vals = ($line =~ /(<$pat(?:\s+.+?)?>.+?<\/$pat>)/g);

      if ( scalar @vals > 0 ) {

        my $curr = 0;
        my $totl = scalar @vals;
        foreach $val (@vals) {
          $curr++;
          ( $tab, $ret ) = process_level ( $val, $pat, $chd, \@group, $tab, $ret, $level - 1, $curr, $totl );
        }

      } else {

        @vals = ($line =~ /(<$pat(?:\s+.+?)\/>)/g);

        my $curr = 0;
        my $totl = scalar @vals;
        foreach $val (@vals) {
          $curr++;
          ( $tab, $ret ) = process_level ( $val, $pat, $chd, \@group, $tab, $ret, $level - 1, $curr, $totl );
        }
      }

    } elsif ( $tmp[$start] eq $capname and $stop - $start >= 2 ) {

      my $pat = $tmp[$start + 1];
      my $chd = "";
      if ( $pat =~ /(.+)\/(.+)/ ) {
        $pat = $1;
        $chd = $2;
      }
      $start += 2;
      my @group = splice (@tmp, $start, $stop - $start);

      # tokenized pattern matching for XML with nested tags

      my @vals = ();
      my $lvl = 0;
      my $dpth = 0;
      my $in_pat = false;
      my @working = ();

      my @tokens = split (/(?<=>)(?=<)/, $line);
      foreach $tkn (@tokens) {

        if ( $lvl < 0 ) {
          # ignore remainder if level drops below zero
        } elsif ( $tkn =~ /^<.+?>.+?<\/.+?>$/ ) {
          # content-containing tag does not change level
          if ( $in_pat ) {
            push (@working, $tkn);
          } elsif ( $tkn =~ /<$pat(?:\s+.+?)?>(.+?)<\/$pat>/ ) {
            # match to full object, push directly to final array
            push (@vals, $tkn);
          }
        } elsif ( $tkn =~ /\/>$/ ) {
          # self-closing token does not change level
          if ( $in_pat ) {
            push (@working, $tkn);
          }
        } elsif ( $tkn =~ /^<[^\/]/ ) {
          # open tag increments level
          if ( $in_pat ) {
            push (@working, $tkn);
          } elsif ( $tkn =~ /<$pat(?:\s+.+?)?>/ ) {
            $in_pat = true;
            $dpth = $lvl;
            push (@working, $tkn);
          }
          $lvl++;
        } elsif ( $tkn =~ /^<\// ) {
          # close tag decrements level
          $lvl--;
          if ( $in_pat ) {
            push (@working, $tkn);
          }
          if ( $tkn =~ /<\/$pat>/ ) {
            if ( $lvl <= $dpth ) {
              $in_pat = false;
              my $str = join ("", @working);
              push (@vals, $str);
              @working = ();
            }
          }
        }
      }
      if ( scalar @working > 0 ) {
        my $str = join ("", @working);
        push (@vals, $str);
        @working = ();
      }

      my $curr = 0;
      my $totl = scalar @vals;
      foreach $val (@vals) {
        $curr++;
        ( $tab, $ret ) = process_level ( $val, $pat, $chd, \@group, $tab, $ret, $level - 1, $curr, $totl );
      }

    } else {

      my @group = splice (@tmp, $start, $stop - $start);
      ( $tab, $ret ) = process_level ( $line, "", "", \@group, $tab, $ret, $level - 1, $indx, $sclr );
    }
  }

  return $tab, $ret;
}

sub process_insd {

  # ... | xtract -insd [+/-] mat_peptide "%peptide" product peptide

  my $quot = shift (@_);
  my @args = @{shift (@_)};

  my $max = scalar @args;

  my @working = ();

  # print common arguments

  push (@working, "-pattern");
  push (@working, "INSDSeq");
  push (@working, "-ACCN");
  push (@working, "INSDSeq_accession-version");

  # collect descriptors

  if ( $args[1] =~ /INSD/ ) {

    push (@working, "-pfx");
    push (@working, "$quot\\n$quot");
    push (@working, "-element");
    push (@working, "$quot\&ACCN$quot");

    for ( my $i = 1; $i < $max; $i++ ) {
      my $val = $args[$i];

      push (@working, "-block");
      push (@working, "INSDSeq");
      push (@working, "-element");
      push (@working, "$quot$val$quot");
    }

    return @working;
  }

  # collect qualifiers

  my @qualifiers = ();
  my @required = ();

  my $feat = "";
  my $partial = 0;

  for ( my $i = 1; $i < $max; $i++ ) {
    my $val = $args[$i];

    if ( $feat eq "" ) {
      if ( $val eq "+" or $val eq "complete" ) {
        $partial = 1;
      } elsif ( $val eq "-" or $val eq "partial" ) {
        $partial = -1;
      } else {
        $feat = $val;
      }
    } else {
      if ( $val =~ /^\+(.+)/ ) {
        $val = $1;
        push (@required, $val);
      }
      push (@qualifiers, $val);
    }
  }

  push (@working, "-group");
  push (@working, "INSDFeature");

  my @itms = split (',', $feat);
  my $fcmd = "-match";
  foreach $itm (@itms) {
    push (@working, "$fcmd");
    push (@working, "INSDFeature_key:$itm");
    $fcmd = "-or";
  }

  if ( $partial == 1 ) {
    push (@working, "-avoid");
    push (@working, "INSDFeature_partial5");
    push (@working, "-and");
    push (@working, "INSDFeature_partial3");
  } elsif ( $partial == -1 ) {
    push (@working, "-match");
    push (@working, "INSDFeature_partial5");
    push (@working, "-or");
    push (@working, "INSDFeature_partial3");
  }

  my $rcmd = "-match";
  foreach my $val (@required) {
    push (@working, "$rcmd");
    push (@working, "INSDQualifier_name:$val");
    $rcmd = "-and";
  }

  push (@working, "-pfx");
  push (@working, "$quot\\n$quot");
  push (@working, "-element");
  push (@working, "$quot&ACCN$quot");

  foreach my $val (@qualifiers) {
    if ( $val =~ /INSD/ ) {
      push (@working, "-block");
      push (@working, "INSDFeature");
      push (@working, "-element");
      push (@working, "$quot$val$quot");

    } else {

      push (@working, "-block");
      push (@working, "INSDQualifier");

      if ( $val =~ /^%(.+)/ ) {
        $val = $1;
        push (@working, "-match");
        push (@working, "INSDQualifier_name:$val");
        push (@working, "-element");
        push (@working, "$quot\%INSDQualifier_value$quot");
      } else {
        push (@working, "-match");
        push (@working, "INSDQualifier_name:$val");
        push (@working, "-element");
        push (@working, "INSDQualifier_value");
      }
    }
  }

  return @working;
}

my $xtract_help = qq{
This Perl version of xtract is obsolete, and is no longer maintained or
supported. The new version, which is written in the Go language, is more
than an order of magnitude faster. To build the new xtract, download and
install the open source Go compiler on your computer, then run:

  go build xtract.go

};

sub xtract {

  # ... | xtract -pattern ... -group ... -block ... -subset ... -element ...

  my $compress_spaces = true;

  my $max = scalar @ARGV;

  if ( $max < 1 ) {
    print STDERR "No arguments were supplied to xtract.pl\n";
    return;
  }

  my $timer = false;

  if ( $ARGV[0] eq "-timer" ) {

    # report execution time of original Perl implementation for comparison to new Go compiled version
    $timer = true;

  } elsif ( $ARGV[0] ne "-fallback" ) {

    # ensure that platform-specific Go compiled version was used if available
    print STDERR "\nPLEASE REWRITE YOUR SCRIPT TO CALL XTRACT INSTEAD OF XTRACT.PL.\n";
    print STDERR "XTRACT WILL RUN AN IMPROVED, COMPILED EXECUTABLE THAT CAN USE\n";
    print STDERR "MULTIPLE CPU CORES TO PROCESS RECORDS BETWEEN ONE AND TWO ORDERS\n";
    print STDERR "OF MAGNITUDE FASTER THAN THE ORIGINAL PERL IMPLEMENTATION.\n\n";
    return;
  }

  # skip past -timer or -fallback command
  shift @ARGV;
  $max = scalar @ARGV;

  if ( $max > 0 and $ARGV[0] eq "-help" ) {
    print "xtract $version (Perl version)\n";
    print $xtract_help;
    return;
  }

  if ( $max > 0 and $ARGV[0] eq "-version" ) {
    print "$version\n";
    return;
  }

  # report deprecated commands

  if ( $max > 0 ) {
    foreach $argument (@ARGV) {
      if ( $argument eq "-present" ) {
        print STDERR "Argument -present is deprecated, use -match Object:Value\n";
      }
      if ( $argument eq "-absent" ) {
        print STDERR "Argument -absent is deprecated, use -avoid Object:Value\n";
      }
      if ( $argument eq "-trim" ) {
        print STDERR "Argument -trim is deprecated, use -block */Object\n";
      }
      if ( $argument eq "-make" or $argument eq "-fuse" or $argument eq "-split" or $argument eq "-pipe" or $argument eq "-repeat" ) {
        print STDERR "Argument '$argument' is deprecated\n";
      }
    }
  }

  # -nocompress allows reading of BLAST XML match data with internal runs of spaces (undocumented)

  if ( $max > 0 and $ARGV[0] eq "-nocompress" ) {
    shift @ARGV;
    $max = scalar @ARGV;
    $compress_spaces = false;
  }

  # -insd to simplify extraction of INSDSeq qualifiers

  if ( $max > 0 and $ARGV[0] eq "-insd" ) {
    if ( $max < 3 and $ARGV[1] !~ /INSD/) {
      print STDERR "Must supply a feature key and at least one qualifier name\n";
      return;

    } elsif ( -t STDIN ) {

      # -insd without piped input will generate an extraction script

      @ARGV = process_insd ( "\"", \@ARGV );
      print "xtract ";
      foreach $itm (@ARGV) {
        print "$itm ";
      }
      print "| \\\n";
      return;

    } else {

      # -insd with piped input will dynamically execute the extraction

      @ARGV = process_insd ( "", \@ARGV );
      $max = scalar @ARGV;
    }
  }

  # recommended nested organizer levels are -group, -block, and -subset

  # intermediate levels (-division, -branch, -section, and -unit) kept in reserve

  # initial capitalized versions (-Group, -Block, etc.) use tokenized pattern matching

  if ( $max > 1 and $ARGV[0] eq "-split" ) {

    # -split reads data stream and writes chunks on separate lines

    my $pat = $ARGV[1];
    my $in_pat = false;
    my @working = ();

    while ( defined($thisline = <STDIN>) ) {
      $thisline =~ s/\r//g;
      $thisline =~ s/\n//g;
      $thisline =~ s/\t//g;
      if ( $compress_spaces ) {
        $thisline =~ s/ +/ /g;
      }
      $thisline =~ s/> +</></g;

      my @tokens = split (/(?<=>)(?=<)/, $thisline);
      foreach $tkn (@tokens) {
        if ( $tkn =~ /^ +(.+)$/ ) {
          $tkn = $1;
        }
        if ( $tkn =~ /^(.+) +$/ ) {
          $tkn = $1;
        }
        if ( $tkn =~ /<$pat(?:\s+.+?)?>(.+?)<\/$pat>/ ) {
          if ( scalar @working > 0 ) {
            my $str = join ("", @working);
            print "$str\n";
            @working = ();
          }
          $in_pat = false;
          print "$tkn\n";
        } elsif ( $tkn =~ /<$pat(?:\s+.+?)?>/ ) {
          push (@working, $tkn);
          $in_pat = true;
        } elsif ( $tkn =~ /<\/$pat>/ ) {
          push (@working, $tkn);
          $in_pat = false;
          my $str = join ("", @working);
          print "$str\n";
          @working = ();
        } elsif ( $in_pat) {
          push (@working, $tkn);
        }
      }
    }
    if ( scalar @working > 0 ) {
      my $str = join ("", @working);
      print "$str\n";
    }

    return;
  }

  if ( $max > 0 and $ARGV[0] eq "-format" ) {

    my $lvl = 0;

    my $xml_ok = true;
    my $doctype_ok = true;

    my $root_token = "";
    my $print_head = false;
    my $print_tail = false;

    my $xml = "";
    my $doctype = "";
    my $pfx = "";
    my $sfx = "";

    my $go_on = true;

    # support -pfx, -sfx, -xml, and -doctype arguments

    for ( my $i = 1; $i < $max; $i++ ) {

      if ( $ARGV[$i] eq "-pfx" ) {
        $i++;
        $pfx = convert_slash ( $ARGV[$i] );
      } elsif ( $ARGV[$i] eq "-sfx" ) {
        $i++;
        $sfx = convert_slash ( $ARGV[$i] );
      } elsif ( $ARGV[$i] eq "-xml" ) {
        $i++;
        $xml = convert_slash ( $ARGV[$i] );
      } elsif ( $ARGV[$i] eq "-doctype" ) {
        $i++;
        $doctype = convert_slash ( $ARGV[$i] );
      }
    }

    while ( $go_on ) {

      $thisline = "";
      if ( $xml ne "" ) {
        $thisline = $xml;
        $xml = "";
      } elsif ( $doctype ne "" ) {
        $thisline = $doctype;
        $doctype = "";
      } elsif ( $pfx ne "" ) {
        $thisline = $pfx;
        $pfx = "";
      } elsif ( defined($thisline = <STDIN>) ) {
      } elsif ( $sfx ne "" ) {
        $thisline = $sfx;
        $sfx = "";
      } else {
        $go_on = false;
      }

      if ( $go_on and $thisline ne "" ) {
        $thisline =~ s/\r//g;
        $thisline =~ s/\n//g;
        $thisline =~ s/\t//g;
        if ( $compress_spaces ) {
          $thisline =~ s/ +/ /g;
        }
        $thisline =~ s/> +</></g;

        my @tokens = split (/(?<=>)(?=<)/, $thisline);
        foreach $tkn (@tokens) {
          if ( $tkn =~ /^ +(.+)$/ ) {
            $tkn = $1;
          }
          if ( $tkn =~ /^(.+) +$/ ) {
            $tkn = $1;
          }
          if ( $lvl < 0 ) {
            # ignore remainder if level drops below zero
          } elsif ( $tkn =~ /^<\?xml/ ) {
            if ( $xml_ok ) {
              # only print first xml line
              print "$tkn\n";
              $xml_ok = false;
            }
          } elsif ( $tkn =~ /^<!DOCTYPE (\S+)/ ) {
            if ( $doctype_ok ) {
              # only print first DOCTYPE line, extract root token
              $root_token = $1;
              print "$tkn\n";
              $doctype_ok = false;
              $print_head = true;
            }
          } elsif ( $tkn =~ /^<\?.+\?>$/ ) {
            # processing instruction
            $lvl++;
            for ( $i = 1; $i < $lvl; $i++ ) {
              print "  ";
            }
            print "$tkn\n";
            $lvl--;
          } elsif ( $tkn =~ /^<\!--.+-->$/ ) {
            # comment
            $lvl++;
            for ( $i = 1; $i < $lvl; $i++ ) {
              print "  ";
            }
            print "$tkn\n";
            $lvl--;
          } elsif ( $tkn =~ /^<.+?>.+?<\/(.+?)>$/ ) {
            # content-containing tag
            if ( $tkn =~ /^<(\S+).*?>.*>$/ ) {
              $lvl++;
              for ( $i = 1; $i < $lvl; $i++ ) {
                print "  ";
              }
              print "$tkn\n";
              $lvl--;
            }
          } elsif ( $tkn =~ /^<.+?\/>$/ ) {
            # self-closing token
            if ( $tkn =~ /^<(\S+).*?\/>$/ ) {
              $lvl++;
              for ( $i = 1; $i < $lvl; $i++ ) {
                print "  ";
              }
              print "$tkn\n";
              $lvl--;
            }
          } elsif ( $tkn =~ /^<[^\/]/ ) {
            # open tag increments level
            $lvl++;
            if ( $tkn eq "<$root_token>" ) {
              if ( $print_head ) {
                print "<$root_token>\n";
                $print_head = false;
                $print_tail = true;
              }
            } else {
              if ( $tkn =~ /^<(\S+).*?>$/ ) {
                for ( $i = 1; $i < $lvl; $i++ ) {
                  print "  ";
                }
                print "$tkn\n";
              } elsif ( $tkn =~ /^<(\S+).*?$/ ) {
                for ( $i = 1; $i < $lvl; $i++ ) {
                  print "  ";
                }
                print "$tkn";
              }
            }
          } elsif ( $tkn =~ /^<\// ) {
            # close tag decrements level
            if ( $tkn ne "</$root_token>" ) {
              for ( $i = 1; $i < $lvl; $i++ ) {
                print "  ";
              }
              print "$tkn\n";
            }
            $lvl--;
          } elsif ( $tkn =~ /^>$/ ) {
            print "$tkn\n";
          } elsif ( $tkn =~ />$/ ) {
            print " $tkn\n";
          } else {
            print " $tkn";
          }
        }
      }
    }

    if ( $print_tail ) {
      print "</$root_token>\n";
    }

    return;
  }

  my $ret = "";
  my $tab = "";

  if ( $max > 2 and $ARGV[0] eq "-pipe" ) {

    # -pipe reads data that has already been compressed by xtract -split

    shift @ARGV;
    while ( defined($thisline = <STDIN>) ) {
      $thisline =~ s/\r//g;
      $thisline =~ s/\n//g;
      $thisline =~ s/\t//g;
      if ( $compress_spaces ) {
        $thisline =~ s/ +/ /g;
      }
      $thisline =~ s/> +</></g;

      my @tmp = @ARGV;
      %memory = ();
      ( $tab, $ret ) = process_level ( $thisline, "", "", \@tmp, "", "", $initial_level, 1, 1 );
      print "$ret";
    }

    return;
  }

  my $begin_time = time();

  if ( $max > 2 and $ARGV[0] eq "-pattern" ) {

    # simple -pattern reads data stream and processes one pattern at a time

    my $pat = $ARGV[1];
    my $chd = "";
    if ( $pat =~ /(.+)\/(.+)/ ) {
      $pat = $1;
      $chd = $2;
    }

    @ARGV = splice (@ARGV, 2, $max - 2);

    my $in_pat = false;
    my @working = ();

    while ( defined($thisline = <STDIN>) ) {
      $thisline =~ s/\r//g;
      $thisline =~ s/\n//g;
      $thisline =~ s/\t//g;
      if ( $compress_spaces ) {
        $thisline =~ s/ +/ /g;
      }
      $thisline =~ s/> +</></g;

      my @tokens = split (/(?<=>)(?=<)/, $thisline);
      foreach $tkn (@tokens) {
        if ( $tkn =~ /^ +(.+)$/ ) {
          $tkn = $1;
        }
        if ( $tkn =~ /^(.+) +$/ ) {
          $tkn = $1;
        }
        if ( $tkn =~ /<$pat(?:\s+.+?)?>/ ) {
          if ( $tkn =~ /<$pat(?:\s+.+?)?>(.+?)<\/$pat>/ ) {
            if ( scalar @working > 0 ) {
              my $str = join ("", @working);
              @working = ();
              my @tmp = @ARGV;
              %memory = ();
              ( $tab, $ret ) = process_level ( $str, $pat, $chd, \@tmp, "", "", $initial_level, 1, 1 );
              print "$ret";
            }
            $in_pat = false;
            my @tmp = @ARGV;
            %memory = ();
            ( $tab, $ret ) = process_level ( $tkn, $pat, $chd, \@tmp, "", "", $initial_level, 1, 1 );
            print "$ret";
          } else {
            push (@working, $tkn);
            $in_pat = true;
          }
        } elsif ( $tkn =~ /<\/$pat>/ ) {
          push (@working, $tkn);
          $in_pat = false;
          my $str = join ("", @working);
          @working = ();
          my @tmp = @ARGV;
          %memory = ();
          ( $tab, $ret ) = process_level ( $str, $pat, $chd, \@tmp, "", "", $initial_level, 1, 1 );
          print "$ret";
        } elsif ( $in_pat) {
          push (@working, $tkn);
        }
      }
    }
    if ( scalar @working > 0 ) {
      my $str = join ("", @working);
      @working = ();
      my @tmp = @ARGV;
      %memory = ();
      ( $tab, $ret ) = process_level ( $str, $pat, $chd, \@tmp, "", "", $initial_level, 1, 1 );
      print "$ret";
    }

    my $end_time = time();
    my $elapsed = $end_time - $begin_time;
    if ( $timer ) {
      if ( $elapsed > 1 ) {
        print STDERR "\nXtract.pl ran in $elapsed seconds\n\n";
      } elsif ( $elapsed > 0 ) {
        print STDERR "\nXtract.pl ran in $elapsed second\n\n";
      }
    }

    return;
  }

  if ( $max > 2 and $ARGV[0] eq "-Pattern" ) {

    # capitalized -Pattern tracks depth of top-level pattern (e.g., for nested Taxon structure)

    my $pat = $ARGV[1];
    my $chd = "";
    if ( $pat =~ /(.+)\/(.+)/ ) {
      $pat = $1;
      $chd = $2;
    }

    @ARGV = splice (@ARGV, 2, $max - 2);

    my $in_pat = false;
    my $lvl = 0;
    my @working = ();

    while ( defined($thisline = <STDIN>) ) {
      $thisline =~ s/\r//g;
      $thisline =~ s/\n//g;
      $thisline =~ s/\t//g;
      if ( $compress_spaces ) {
        $thisline =~ s/ +/ /g;
      }
      $thisline =~ s/> +</></g;

      my @tokens = split (/(?<=>)(?=<)/, $thisline);
      foreach $tkn (@tokens) {
        if ( $tkn =~ /^ +(.+)$/ ) {
          $tkn = $1;
        }
        if ( $tkn =~ /^(.+) +$/ ) {
          $tkn = $1;
        }
        if ( $tkn =~ /<$pat(?:\s+.+?)?>/ ) {
          if ( $tkn =~ /<$pat(?:\s+.+?)?>(.+?)<\/$pat>/ ) {
            if ( scalar @working > 0 ) {
              my $str = join ("", @working);
              @working = ();
              my @tmp = @ARGV;
              %memory = ();
              ( $tab, $ret ) = process_level ( $str, $pat, $chd, \@tmp, "", "", $initial_level, 1, 1 );
              print "$ret";
            }
            $in_pat = false;
            my @tmp = @ARGV;
            %memory = ();
            ( $tab, $ret ) = process_level ( $tkn, $pat, $chd, \@tmp, "", "", $initial_level, 1, 1 );
            print "$ret";
          } else {
            $lvl++;
            push (@working, $tkn);
            $in_pat = true;
          }
        } elsif ( $tkn =~ /<\/$pat>/ ) {
          $lvl--;
          if ($lvl < 1) {
            push (@working, $tkn);
            $in_pat = false;
            my $str = join ("", @working);
            @working = ();
            my @tmp = @ARGV;
            %memory = ();
            ( $tab, $ret ) = process_level ( $str, $pat, $chd, \@tmp, "", "", $initial_level, 1, 1 );
            print "$ret";
          } elsif ( $in_pat) {
            push (@working, $tkn);
          }
        } elsif ( $in_pat) {
          push (@working, $tkn);
        }
      }
    }
    if ( scalar @working > 0 ) {
      my $str = join ("", @working);
      @working = ();
      my @tmp = @ARGV;
      %memory = ();
      ( $tab, $ret ) = process_level ( $str, $pat, $chd, \@tmp, "", "", $initial_level, 1, 1 );
      print "$ret";
    }

    return;
  }

  my $rpt = 1;

  # read entire XML input stream into a single string

  $holdTerminator = $/;
  undef $/;
  $data = <STDIN>;
  $/ = $holdTerminator;

  # remove newlines, tabs, space between tokens, compress runs of spaces,

  $data =~ s/\r//g;
  $data =~ s/\n//g;
  $data =~ s/\t//g;
  if ( $compress_spaces ) {
   $data =~ s/ +/ /g;
  }
  $data =~ s/> +</></g;

  if ( $max == 0 or ( $max == 1 and $ARGV[0] ne "-outline" and $ARGV[0] ne "-synopsis" )) {

    # no useful arguments, just print entire line of compressed XML

    print "$data\n";

    return;
  }

  if ( $max > 2 and $ARGV[0] eq "-repeat" ) {

    # first look for -repeat count for performance testing

    $rpt = $ARGV[1];
    @ARGV = splice (@ARGV, 2, $max - 2);
    $max = scalar @ARGV;
  }

  if ( $ARGV[0] ne "-pattern" ) {

    # no -pattern, process entire XML without looking for pattern

    my @args = @ARGV;
    %memory = ();
    ( $tab, $ret ) = process_level ( $data, "", "", \@args, "", "", $initial_level, 1, 1 );
    print "$ret";

    return;
  }

  # split by -pattern - can handle single term or parent/child (including heterogeneous parent/*)

  my $ptrn = $ARGV[1];
  my $chld = "";
  if ( $ptrn =~ /(.+)\/(.+)/ ) {
    $ptrn = $1;
    $chld = $2;
  }
  @vals = ($data =~ /(<$ptrn(?:\s+.+?)?>.+?<\/$ptrn>)/g);

  if ( $max == 2 ) {
    foreach $val (@vals) {
      print "$val\n";
    }

    return;
  }


  # process remaining arguments

  my @args = splice (@ARGV, 2, $max - 2);

  my $curr = 0;
  my $totl = scalar @vals;

  for ( my $idx = 0; $idx < $rpt; $idx++ ) {
    foreach $val (@vals) {
      %memory = ();
      ( $tab, $ret ) = process_level ( $val, $ptrn, $chld, \@args, "", "", $initial_level, $curr, $totl );
      print "$ret";
    }
  }
}

# execute XML extraction

xtract ();

# if -synopsis, print cumulative summary of XML paths

# if ( scalar (keys %synopsis_acc) > 0 ) {
#   my @keys;
#   @keys = keys %synopsis_acc;
#   for my $ky (sort @keys) {
#     $vl = $synopsis_acc{$ky};
#     for ( $i = 0; $i < $vl; $i++ ) {
#     	print "$ky\n";
#     }
#   }
# }

# close input and output files

close (STDIN);
close (STDOUT);
close (STDERR);
