
import java.awt.BasicStroke;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import javax.imageio.ImageIO;

public final class Tandems {

    public void SetSequences(String[] seq, String[] sname) {
        this.seq = seq;
        this.sname = sname;
        nseq = seq.length;
    }

    public void SetFileName(String a) {
        filePath = a;
    }

    public void SetKmerLen(int i) {
        kmerln = i;
        if (kmerln < 12) {
            kmerln = 12;
        }
        minlenblock = kmerln + kmerln;
        minlenseq = minlenblock;
    }

    public void SetImage(int w, int h) {
        if (w > 0) {
            iwidth = w;
        }
        if (h > 0) {
            iheight = h;
        }
    }

    public void SetNumberUnits(int i) {
        nblocks = i - 1;
        if (nblocks < 1) {
            nblocks = 1;
        }
        if (nblocks > 1000) {
            nblocks = 1000;
        }
    }

    public void SetTandemLen(int i) {
        minlenblock = i;
        if (minlenblock < mnblock) {
            minlenblock = mnblock;
        }
        minlenseq = minlenblock;
    }

    public void SetSequenceLen(int i) {
        minlenseq = i;
        if (minlenseq < minlenblock) {
            minlenseq = minlenblock;
        }
    }

    public void SetTelomers(boolean i) {
        if (i) {
            telomers = 14;
        }
    }

    public void SetFlanks(int i) {
        flanks = i;
    }

    public void SetShowSeq(boolean i) {
        SeqShow = i;
    }

    public void Run() throws IOException {
        startTime = System.nanoTime();
        for (int i = 0; i < nseq; i++) {
            ArrayList<Integer> u = Mask(seq[i], kmerln, kmerln);
            ClusteringMasking(seq[i], u, kmerln);
            PrintResult(i);
            PictureSave(i, iwidth, iheight);
        }
    }

    public int[] ArrayExtend(int[] srcArray, int n) {
        int[] destArray = new int[srcArray.length + n];
        System.arraycopy(srcArray, 0, destArray, 0, srcArray.length);
        return destArray;
    }

    public int[] ArrayTrim(int[] srcArray, int n) {
        int[] destArray = new int[n];
        System.arraycopy(srcArray, 0, destArray, 0, n);
        return destArray;
    }

    public byte[] ArrayExtendByte(byte[] srcArray, int n) {
        byte[] destArray = new byte[srcArray.length + n];
        System.arraycopy(srcArray, 0, destArray, 0, srcArray.length);
        return destArray;
    }

    private ArrayList<Integer> Mask(String seq, int kmer, int kmer2) {
        int l = seq.length();
        int i = 0;
        int j = 0;
        int p = 0;

        String s;
        String aseq = dna.ComplementDNA(seq);
        HashMap<String, Integer> px2 = new HashMap<>();

        int[] ax = new int[kmer];
        int[] bx = new int[5];

        byte b[] = seq.getBytes();
        for (i = 0; i < l; i++) {
            b[i] = tables.dx2[b[i]];
        }
        b = ArrayExtendByte(b, l + 1);
        b[l] = 4;
        for (i = 1; i < l + 1; i++) {
            b[l + i] = tables.cdnat2[b[l - i]];
        }
        for (i = 0; i < kmer - 1; i++) {
            ax[i] = b[i];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[i];
            bx[ax[kmer - 1]]++;
            if (bx[4] == 0) {
                s = seq.substring(i - kmer + 1, i + 1);
                if (px2.containsKey(s)) {
                    p = px2.get(s) + 1;
                    px2.put(s, p);
                } else {
                    px2.put(s, 1);
                }
            }
            bx[b[i + 1 - kmer]]--;
            for (j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }
        //reverse
        bx = new int[5];
        for (i = 0; i < kmer - 1; i++) {
            ax[i] = b[l + i + 1];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[l + i + 1];
            bx[ax[kmer - 1]]++;
            if (bx[4] == 0) {
                s = aseq.substring(i - kmer + 1, i + 1);
                if (px2.containsKey(s)) {
                    p = px2.get(s) + 1;
                    px2.put(s, p);
                }
            }
            bx[b[l + i + 2 - kmer]]--;
            for (j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }

        int[] u = new int[l];
        bx = new int[5];
        for (i = 0; i < kmer - 1; i++) {
            ax[i] = b[i];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[i];
            bx[ax[kmer - 1]]++;
            int x = i - kmer + 1;
            if (bx[4] == 0) {
                s = seq.substring(x, i + 1);
                p = px2.get(s);
                if (p > nblocks) {
                    for (j = x; j < x + kmer; j++) {
                        u[j] = 1;
                    }
                }
            }
            bx[b[i + 1 - kmer]]--;
            for (j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }
        //reverse
        bx = new int[5];
        for (i = 0; i < kmer - 1; i++) {
            ax[i] = b[l + i + 1];
            bx[ax[i]]++;
        }
        for (i = kmer - 1; i < l; i++) {
            ax[kmer - 1] = b[l + i + 1];
            bx[ax[kmer - 1]]++;
            if (bx[4] == 0) {
                int x = i - kmer + 1;
                s = aseq.substring(x, i + 1);
                if (px2.containsKey(s)) {
                    p = px2.get(s);
                    if (p > nblocks) {
                        for (j = x; j < x + kmer; j++) {
                            u[j] = 2;
                        }
                    }
                }
            }
            bx[b[l + i + 2 - kmer]]--;
            for (j = 0; j < kmer - 1; j++) {
                ax[j] = ax[j + 1];
            }
        }

        ArrayList<Integer> z = new ArrayList<>();
        for (i = 0; i < l; i++) {
            if (u[i] > 0) {
                int x1 = i;
                int x2 = i;
                for (j = i + 1; j < l; j++) {
                    if (u[j] > 0) {
                        x2 = j;
                    } else {
                        i = j;
                        z.add(x1);
                        z.add(x2);
                        break;
                    }
                }
            }
        }
        for (i = 0; i < z.size(); i += 2) {
            int x1 = z.get(i);
            int x2 = z.get(i + 1);

            if (x1 > -1) {
                int w = x2;
                for (j = i + 2; j < z.size(); j += 2) {
                    int x3 = z.get(j);
                    int x4 = z.get(j + 1);
                    if (x3 > -1) {
                        if (w + kmer2 > x3) {
                            w = x4;
                            x2 = x4;
                            z.set(i + 1, x4);
                            z.set(j, -1);
                            z.set(j + 1, -1);
                        }
                    } else {
                        break;
                    }
                }
            }
        }

        ArrayList<Integer> z2 = new ArrayList<>();

        for (i = 0; i < z.size(); i += 2) {
            int x1 = z.get(i);
            if (x1 > -1) {
                int x2 = z.get(i + 1);
                if (x2 - x1 > kmer2) {
                    z2.add(x1);
                    z2.add(x2);
                }
            }
        }
        return z2;
    }

    private int ClusteringMasking(String seq, ArrayList<Integer> u, int kmer) {
        int n = u.size();
        if (n < 2) {
            return -1;
        }
        ArrayList<String> z = new ArrayList<>();
        ArrayList<Integer> x = new ArrayList<>();
        for (int i = 0; i < n; i = i + 2) {
            int x1 = u.get(i);
            int x2 = u.get(i + 1);
            String s = seq.substring(x1, 1 + x2);
            x.add(x1);// location on sequence
            z.add(s); // repeat sequence
        }

        String[] s = new String[z.size()];
        s = z.toArray(s);

        SequencesClustering sc = new SequencesClustering(s, 90, true, kmer, minlenseq);
        int[] q = sc.Result();
        int ncl = sc.getNcl();

        if (ncl < 1) {
            return -1;
        }

        bb = new ArrayList<>();
        for (int f = 1; f < ncl + 1; f++) {
            int[] k7 = new int[ncl + ncl + 1];
            int h = 0;
            k7[0] = 0;
            for (int j = 0; j < q.length; j++) {
                if (q[j] == f) {
                    h++;
                    k7[h] = x.get(j);
                    h++;
                    k7[h] = s[j].length();
                    k7[0] = k7[0] + 2;
                }
            }
            bb.add(ArrayTrim(k7, h + 1));
        }
        return bb.size();
    }

    private void PrintResult(int n) throws IOException {
        if (bb == null) {
            return;
        }
        int k = 0;
        double z = numnonn;
        int l = seq[n].length();
        int v = numnonn;
        long duration = (System.nanoTime() - startTime) / 1000000000;

        String reportfile = filePath + "_" + (n + 1) + ".gff";
        String maskedfile = filePath + "_" + (n + 1) + ".msk";
        if (nseq == 1) {
            reportfile = filePath + ".gff";
            maskedfile = filePath + ".msk";
        }

        byte[] m = new byte[l];

        try (FileWriter fileWriter = new FileWriter(maskedfile)) {
            fileWriter.write(">" + sname[n]);
            System.out.println("Saving masked file: " + maskedfile);
            for (int i = 0; i < bb.size(); i++) {
                int[] z7 = bb.get(i);
                if (z7[0] > 1) {
                    for (int j = 1; j < z7[0]; j += 2) {
                        int x = z7[j] + Math.abs(z7[j + 1]) - 1;
                        for (int h = z7[j]; h < x; h++) {
                            m[h] = 4;
                        }
                    }
                }
            }
            char[] c = seq[n].toCharArray();
            for (int i = 0; i < l; i++) {
                if (m[i] == 0) {
                    z--;
                    c[i] = (char) (c[i] - 32);
                }
            }
            z = (z * 100 / v);
            fileWriter.write("Sequence coverage by repeats " + String.format("%.2f", z) + "%\n\n");
            fileWriter.write(new String(c));
            System.out.println("Sequence coverage by repeats " + String.format("%.2f", z) + "%");
        }

        StringBuilder sr = new StringBuilder();
        sr.append("kmer=").append(kmerln).append("\n").append("Initial string length filter=").append(minlenblock).append("\n").append("String length filter=").append(minlenseq).append("\n").append("\n");
        sr.append("__________________________________________________\n Repeats search for: ").append(filePath).append("//").append(sname[n]).append(" ").append(l).append("bp :\n");
        sr.append("Time taken: ").append(duration).append(" seconds\n\n");
        if (SeqShow) {
            sr.append("Seqid\tRepeat\tClusterID\tStart\tStop\tLength\tStrand\tPhase\tSequence").append("\n");
        } else {
            sr.append("Seqid\tRepeat\tClusterID\tStart\tStop\tLength\tStrand\tPhase\tSequence").append("\n");
        }
        if (v < l) {
            double d = ((l - v) * 100) / l;
            sr.append("Sequence gap (bp)=").append(l - v).append(" (").append(String.format("%.3f", d)).append("%)\n");
        }

        try (FileWriter fileWriter = new FileWriter(reportfile)) {
            System.out.println("Saving report file: " + reportfile);
            fileWriter.write("Sequence coverage by repeats " + String.format("%.2f", z) + "%\n\n");

            fileWriter.write(sr.toString());
            for (int i = 0; i < bb.size(); i++) {
                int[] z7 = bb.get(i);
                if (z7[0] > 1) {
                    k++;
                    for (int j = 1; j < z7[0]; j += 2) {
                        if (Math.abs(z7[j + 1]) > minlenseq) {

                            String s0 = "";
                            int x = z7[j] + Math.abs(z7[j + 1]) - 1;

                            if (SeqShow) {
                                if (x > seq[n].length()) {
                                    s0 = seq[n].substring(z7[j]);
                                } else {
                                    s0 = seq[n].substring(z7[j], x);
                                }
                                if (flanks > 0) {
                                    String s1 = "";
                                    String s2 = "";
                                    if (z7[j] - flanks > 0) {
                                        s1 = seq[n].substring(z7[j] - flanks, z7[j]).toUpperCase();
                                    } else {
                                        if (z7[j] > 1) {
                                            s1 = seq[n].substring(1, z7[1] - 1).toUpperCase();
                                        }
                                    }
                                    if (x + flanks < l) {
                                        s2 = seq[n].substring(x, x + flanks).toUpperCase();
                                    } else {
                                        if (l - x > 0) {
                                            s2 = seq[n].substring(x, l).toUpperCase();
                                        }
                                    }
                                    s0 = s1 + s0 + s2;
                                }
                                if (z7[j + 1] < 0) {
                                    s0 = dna.ComplementDNA2(s0);
                                }
                            }
                            sr = new StringBuilder();
                            if (z7[j + 1] > 0) {
                                sr.append(sname[n]).append("\t").append(".").append("\t").append(k).append("\t").append(z7[j] + 1).append("\t").append(x + 1).append("\t").append(z7[j + 1]).append("\t").append(".").append("\t").append("+").append("\t").append(s0).append("\n");
//                                sr.append(k).append("\t").append(z7[j] + 1).append("\t").append(x + 1).append("\t").append(z7[j + 1]).append("\t").append("\t").append(s0).append("\n");
                                fileWriter.write(sr.toString());
                            } else {
                                sr.append(sname[n]).append("\t").append(".").append("\t").append(k).append("\t").append(-z7[j + 1]).append("\t").append(x + 1).append("\t").append(z7[j] - 1).append("\t").append(".").append("\t").append("-").append("\t").append(s0).append("\n");
//                               sr.append(k).append("\t").append(x + 1).append("\t").append(-z7[j] - 1).append("\t").append(-z7[j + 1]).append("\t").append(s0).append("\n");
                                fileWriter.write(sr.toString());
                            }
                        }
                    }
                }
            }
        }
        /*
GFF format General Feature Format is a format for describing genes and other features associated with DNA, RNA and Protein sequences. GFF lines have nine tab-separated fields:
Generic Feature Format Version 3 (GFF3) https://github.com/The-Sequence-Ontology/Specifications/blob/master/gff3.md
1. seqid - Must be a chromosome or scaffold or contig.
2. source - The program that generated this feature.
3. type - "repeat".
4. start - The starting position of the feature in the sequence. The first base is numbered 1.
5. stop - The ending position of the feature (inclusive).
6. score - length 
7. strand - Valid entries include '+', '-', or '.' (for don't know/care).
8. phase - If the feature is a coding exon, frame should be a number between 0-2 that represents the reading frame of the first base. If the feature is not a coding exon, the value should be '.'.
9. attributes - All lines with the same group are linked together into a single item.
         */
    }

    private void PictureSave(int n, int dw, int dh) throws IOException {
        if (bb == null) {
            return;
        }
        int k = 100;
        int b = bb.size();
        int l = seq[n].length();
        int width = l < 5_000_000 ? l / 100 : 5_000 + (l - 5_000) / 200; // image=10000x3000
        if (dw > 0) {
            width = dw;
        }
        if (width > 46340) { //too big a picture leads to an error  
            width = 46340; //https://jobcardsystems.com/index.php/blog/29-effective-handling-of-big-images-in-java
        }
        if (width < 500) {
            width = 500;
        }
        int height = b < 5000 ? k + b : 5100 + (b - 5000) / 200;        //too big a picture leads to an error
        if (dh > 0) {
            height = dh;
        }
        if (height > 46340) {
            height = 46340;
        }
        if (height < 500) {
            height = 500;
        }
        float dotSize = 1 + (b / 500);                                  //7.0f;
        double w1 = (double) width / l;                                 // nucleotides per pixel        
        if (dotSize < 1) {
            dotSize = 1.0f;
        }
        if (dotSize > 10) {
            dotSize = 10.0f;
        }

        String pngfile = filePath + "_" + (n + 1) + ".png";
        if (nseq == 1) {
            pngfile = filePath + ".png";
        }

        BufferedImage image = new BufferedImage(width + 100, height + 100, BufferedImage.TYPE_INT_RGB);
        Graphics2D g2d = image.createGraphics();

        g2d.setStroke(new BasicStroke(dotSize));
        g2d.setColor(Color.WHITE);
        g2d.fillRect(0, 0, width + 100, height + 100);
        g2d.setColor(Color.BLACK);
        g2d.drawLine(50, 15, width + 50, 15);

        int w = (width + 100) / 10;
        int d = l / 10;
        for (int i = 0; i < 11; i++) {
            g2d.drawLine(i * w + 50, 5, i * w + 50, 20);
            g2d.drawString(String.format("%,d", (1 + (int) (i * d))), 55 + i * w, 10);
        }

        for (int i = 0; i < b; i++) {
            int[] z7 = bb.get(i);
            if (z7[0] > 1) {
                int y = (i + 70);// y1-y2 line                
                for (int j = 1; j < z7[0]; j += 2) {
                    if (Math.abs(z7[j + 1]) > minlenseq) {
                        int x1 = 50 + (int) (z7[j] * w1);
                        int x2 = 50 + (int) ((z7[j] + Math.abs(z7[j + 1])) * w1);

                        if (z7[j + 1] > 0) {
                            g2d.setColor(Color.BLUE);
                        } else {
                            g2d.setColor(Color.RED);
                        }
                        g2d.drawLine(x1, y, x2, y); //(int x1, int y1, int x2, int y2)
                    }
                }
            }
        }
        g2d.dispose();
        try {
            System.out.println("Saving picture : " + pngfile);
            File outputFile = new File(pngfile);
            ImageIO.write(image, "png", outputFile);
        } catch (IOException e) {
            System.out.println("Saving picture is large.");
        }
    }

    private String filePath;
    public int nseq;
    public int iwidth = 0;
    public int iheight = 0;
    private String[] seq;
    private String[] sname;
    private int nblocks = 1;
    private int numnonn = 0;        // calculated non-N bases at the sequence
    private final int ssrlen = 30;
    private int minlenblock = 50;   // initial sequence length user control  
    private int minlenseq = 50;     // sequence length 
    private final int mnblock = 20;
    private int kmerln = 21;
    private int flanks = 20;
    private int telomers = 11;      // Kmax=11 ->SSR  Kmax=14 -> telomers
    private boolean SeqShow;
    private ArrayList<int[]> bb;
    private long startTime;
}
