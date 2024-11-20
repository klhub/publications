a=[ 1 2 3 4]';
a=a/norm(a)
a1=a(1);a2=a(2);a3=a(3);a4=a(4);
b=[ 1 2 3 4]';
b=b/norm(b)
b1=b(1);b2=b(2);b3=b(3);b4=b(4);

a = [ +a4 +a3 -a2 +a1
      -a3 +a4 +a1 +a2
      +a2 -a1 +a4 +a3
      -a1 -a2 -a3 +a4 ] ;
a*b