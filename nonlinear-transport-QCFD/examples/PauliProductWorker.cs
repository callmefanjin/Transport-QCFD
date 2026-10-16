// Exact ordered single-Pauli rotations, evaluated as in-place amplitude pairs.
// Runtime helper only; no alteration of the solver's splitting or coefficients.
using System;
using System.IO;
class PauliProductWorker {
    static int Parity(int x) { x ^= x >> 16; x ^= x >> 8; x ^= x >> 4; return (0x6996 >> (x & 15)) & 1; }
    static void Main(string[] args) {
        int n, m, repetitions;
        int[] xs, zs; double[] cs, fr, fi, re, im;
        using(var b = new BinaryReader(File.OpenRead(args[0]))) {
            n=b.ReadInt32(); m=b.ReadInt32(); repetitions=b.ReadInt32();
            xs=new int[m]; zs=new int[m]; cs=new double[m]; fr=new double[m]; fi=new double[m];
            for(int k=0;k<m;k++) xs[k]=b.ReadInt32();
            for(int k=0;k<m;k++) zs[k]=b.ReadInt32();
            for(int k=0;k<m;k++) cs[k]=b.ReadDouble();
            for(int k=0;k<m;k++) fr[k]=b.ReadDouble();
            for(int k=0;k<m;k++) fi[k]=b.ReadDouble();
            re=new double[n]; im=new double[n];
            for(int i=0;i<n;i++) re[i]=b.ReadDouble();
            for(int i=0;i<n;i++) im[i]=b.ReadDouble();
        }
        var left=new int[n][]; var right=new int[n][]; var signs=new sbyte[n][];
        for(int k=0;k<m;k++) {
            int x=xs[k], z=zs[k];
            if(x!=0 && left[x]==null) {
                left[x]=new int[n/2]; right[x]=new int[n/2]; int p=0;
                for(int i=0;i<n;i++) if(i<(i^x)) {left[x][p]=i;right[x][p++]=i^x;}
            }
            if(signs[z]==null) { signs[z]=new sbyte[n]; for(int i=0;i<n;i++) signs[z][i]=(sbyte)(1-2*Parity(z&i)); }
        }
        for(int step=0;step<repetitions;step++) for(int q=0;q<2*m-1;q++) {
            int k=q<m?q:2*m-2-q; int x=xs[k]; var s=signs[zs[k]]; double c=cs[k];
            if(x==0) {
                for(int i=0;i<n;i++) {double a=re[i],b=im[i],f=fi[k]*s[i];re[i]=c*a-f*b;im[i]=c*b+f*a;}
            } else {
                var l=left[x];var h=right[x];
                if(fi[k]==0) {
                    double factor=fr[k];
                    for(int p=0;p<l.Length;p++) {
                        int i=l[p],j=h[p];double f=factor*s[i],a=re[i],b=im[i],d=re[j],e=im[j];
                        re[i]=c*a+f*d;im[i]=c*b+f*e;re[j]=c*d-f*a;im[j]=c*e-f*b;
                    }
                } else {
                    double factor=fi[k];
                    for(int p=0;p<l.Length;p++) {
                        int i=l[p],j=h[p];double f=factor*s[i],a=re[i],b=im[i],d=re[j],e=im[j];
                        re[i]=c*a-f*e;im[i]=c*b+f*d;re[j]=c*d-f*b;im[j]=c*e+f*a;
                    }
                }
            }
        }
        using(var b=new BinaryWriter(File.Create(args[1]))) {for(int i=0;i<n;i++)b.Write(re[i]);for(int i=0;i<n;i++)b.Write(im[i]);}
    }
}
