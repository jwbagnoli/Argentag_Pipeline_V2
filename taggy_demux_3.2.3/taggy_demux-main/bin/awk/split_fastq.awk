BEGIN {
	reverse["A"] = "T";
	reverse["C"] = "G";
	reverse["G"] = "C";
	reverse["T"] = "A";
	split(brk_offset, offsets, ",");
}

NR==FNR {
	#printf("Original col 2: %s.\n",$2);
	sub(/CS1_EPF_F_/,"",$2);
	switch($4){
	case "fwd":
	  offset = offsets[$2] + 0;
	  break
	case "rev":
	  offset = -offsets[$2] + 0;
	  break
	default:
	  printf("Unexpected orientation %s. Expected \"fwd\" or \"rev\"\n",$4) > "/dev/stderr";
	  exit -1;
	  break;
	}
	#printf("New col 2: %s. Offset: %d.\n",$2,offset);
	brk = $3 + offset + 0;
	b["@" $1] = (b["@" $1])?(b["@" $1] " " brk):brk;
	next
}

(NR!=FNR) && (FNR%4 != 1){
	next
}

(NR!=FNR){
	split($0, a, " ");
	getline seq;
	getline sep;
	getline qual;
	if(a[1] in b){
	  #print b[a[1]]
	  split(b[a[1]], c, " ");
	  n_hits=asort(c,d);
	  right=length(seq);
	  for (hit = n_hits; hit >= 0; hit--){
	    left = (hit == 0)? 1 : d[hit] + 1;
	    #printf("%s: %d-%d\n", hit, left, right);
	    sequence = substr(seq,left,right-left+1)
	    quality = substr(qual,left,right-left+1)
	    printf("%s-%d\n%s\n%s\n%s\n",a[1],hit+1,sequence,sep,quality)
	    right=left-1;
	  }
	}else{
	  printf("%s\n%s\n%s\n%s\n",a[1],seq,sep,qual)
	}
}
