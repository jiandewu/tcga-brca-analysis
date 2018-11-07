#!/usr/bin/awk -f

{
    z1=0;z2=0
    for (i=2;i<=111;i++) 
        if ($i=="0.0") z1++
    if (z1>70)
        for (i=112;i<=NF;i++)
            if ($i=="0.0") z2++
    if (z2<70)
        print
}