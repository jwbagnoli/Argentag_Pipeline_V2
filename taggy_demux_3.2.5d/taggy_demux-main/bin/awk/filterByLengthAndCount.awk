BEGIN{
  if (! min) min = 0;
}

FNR%4 == 1 {
  getline seq;
  getline dummy;
  getline qual;
  len = length(seq);
  if ((len >= min) && ((len <= max) || (!max))) {
    printf("%s\n%s\n+\n%s\n",$0,seq,qual);
    lines++;
  }
}

END {
  printf("%d\n", lines) >"/dev/fd/3"
}
