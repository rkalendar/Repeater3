
import java.io.IOException;
import java.io.File;
import java.nio.file.Files;
import java.nio.file.Paths;

public class Repeater {

    public static void main(String[] args) {
        if (args.length > 0) {
            String infile = args[0]; // file path or Folder
            String s = String.join(" ", args).toLowerCase() + " ";
            int kmer = 21;
            int minlen = 50;
            int seqlen = 100;
            int width = 0;
            int hight = 0;
            int flanksshow = 0;
            boolean seqshow = false;
            boolean ssrrun = false;

            System.out.println("Current Directory: " + System.getProperty("user.dir"));
            System.out.println("Command-line arguments:");
            System.out.println("Target file or Folder: " + infile);

            if (s.contains("seqshow=true")) {
                seqshow = true;
            }

            if (s.contains("flanks=")) {
                flanksshow = 0;
                int j = s.indexOf("flanks=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    flanksshow = StrToInt(s.substring(j + 7, x));
                }
                if (flanksshow < 0) {
                    flanksshow = 0;
                }
                if (flanksshow > 1000) {
                    flanksshow = 1000;
                }
            }
            if (s.contains("kmer=")) {
                int j = s.indexOf("kmer=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    kmer = StrToInt(s.substring(j + 5, x));
                }
                if (kmer < 12) {
                    kmer = 12;
                }
            }
            kmer = 21;
            if (s.contains("image=")) { // image=10000x3000
                int j = s.indexOf("image=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    String[] d = s.substring(j + 6, x).split("x");
                    if (d.length > 0) {
                        width = StrToInt(d[0]);
                        hight = StrToInt(d[1]);
                        if (width < 1000) {
                            width = 1000;
                        }
                        if (hight < 500) {
                            hight = 500;
                        }
                        if (hight == 0 | width == 0) {
                            width = 0;
                            hight = 0;
                        }
                    }
                }
            }

            if (s.contains("min=")) {
                int j = s.indexOf("min=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    minlen = StrToInt(s.substring(j + 4, x));
                }
                if (minlen < kmer) {
                    minlen = kmer;
                }
            }

            if (s.contains("sln=")) {
                int j = s.indexOf("sln=");
                int x = s.indexOf(" ", j);
                if (x > j) {
                    seqlen = StrToInt(s.substring(j + 4, x));
                }
                if (seqlen < minlen) {
                    seqlen = minlen;
                }
            }

            File folder = new File(infile);
            if (folder.exists() && (folder.isDirectory() || folder.isFile())) {
                if (folder.isDirectory()) {
                    File[] files = folder.listFiles();
                    int k = -1;
                    String[] filelist = new String[files.length];
                    for (File file : files) {
                        if (file.isFile()) {
                            filelist[++k] = file.getAbsolutePath();
                        }
                    }
                    for (String nfile : filelist) {
                        try {
                            SaveResult(nfile, kmer, minlen, seqlen, flanksshow, seqshow, ssrrun, width, hight);
                        } catch (Exception e) {
                            System.err.println("Failed to open file: " + nfile);
                        }
                    }

                } else {
                    SaveResult(infile, kmer, minlen, seqlen, flanksshow, seqshow, ssrrun, width, hight);
                }
            }

        } else {
            System.out.println("REPEATER 3 (2024) by Ruslan Kalendar (ruslan.kalendar@helsinki.fi)\nhttps://github.com/rkalendar/Repeater\n");
            System.out.println("Basic usage:");
            System.out.println("java -jar \\Repeater3\\dist\\Repeater3.jar <inputfile>/<inputfolderpath> <optional_commands>");
            System.out.println("Common options:");
            System.out.println("kmer=12\tminimal kmer=12 (default 21)");
            System.out.println("min=100\tinitial repeat length (default min=50), it can be equal to kmer=");
            System.out.println("sln=150\tstring length (default sln=100), it can be equal to min=");
            System.out.println("flangs=100\textend the flanks of the repeat with an appropriate length (100 nt) (default flangs=0)");
            System.out.println("image=10000x3000\t (by default, the dimensionality of the image is automatically determined)");
            System.out.println("mask=true/false\tgenerate a new file with masking repeats (default mask=true)");
            System.out.println("seqshow=true/false\textract repeat sequences (default seqshow=false)");
            System.out.println("quick=true/false\tquick analysis of repeats, without deep analysis and their clustering (default quick=true)");
            System.out.println("ssr=true\tanalyzing only the SSR/telomers loci (default ssr=false)\n");
            System.out.println("java -jar \\Repeater3\\dist\\Repeater3.jar <inputfile> ssr=true seqshow=true flanks=100");
            System.out.println("java -jar \\Repeater3\\dist\\Repeater3.jar <inputfile> min=100 sln=300 quick=false mask=true\n");
            System.out.println("java -jar \\Repeater3\\dist\\Repeater3.jar <inputfile> ssr=true seqshow=true flanks=100");
            System.out.println("Large genome settings:");
            System.out.println("java -jar -Xms16g -Xmx32g \\Repeater3\\dist\\Repeater3.jar <inputfile> kmer=21 quick=false mask=true\n");
            System.out.println("Analysing all files in the directory:");
            System.out.println("java -jar \\Repeater3\\dist\\Repeater3.jar \\Repeater2\\test\\ kmer=21 image=10000x3000 quick=false\n");
        }
    }

    public static int StrToInt(String str) {
        StringBuilder r = new StringBuilder();
        int z = 0;
        r.append(0);
        for (int i = 0; i < str.length(); i++) {
            char chr = str.charAt(i);
            if (chr > 47 && chr < 58) {
                r.append(chr);
                z++;
                if (z > 10) {
                    break;
                }
            }
            if (chr == '.' || chr == ',') {
                break;
            }
        }
        return (Integer.parseInt(r.toString()));
    }

    private static void SaveResult(String infile, int kmer, int minlen, int seqlen, int flanksshow, boolean seqshow, boolean ssrrun, int width, int hight) {
        try {
            long startTime = System.nanoTime();
            byte[] binaryArray = Files.readAllBytes(Paths.get(infile));
            ReadingSequencesFiles rf = new ReadingSequencesFiles(binaryArray);
            if (rf.getNseq() == 0) {
                System.out.println("There is no sequence(s).");
                System.out.println("File format in Fasta:\n>header\nsequence here\n\nIn FASTA format the line before the nucleotide sequence, called the FASTA definition line, must begin with a carat (\">\"), followed by a unique SeqID (sequence identifier).\nThe line after the FASTA definition line begins the nucleotide sequence.\n");
                System.out.println(">seq1\nactacatactacatcactctctctccgcacag\n");
                return;
            }
            System.out.println("Running...");
            Tandems s2 = new Tandems();
            s2.SetSequences(rf.getSequences(), rf.getNames());
            s2.SetKmerLen(kmer);
            s2.SetTandemLen(minlen);
            s2.SetSequenceLen(seqlen);
            s2.SetShowSeq(seqshow);
            s2.SetFlanks(flanksshow);
            s2.SetFileName(infile);
            System.out.println("Target file: " + infile);
            System.out.println("Target sequence length = " + rf.getLength() + " nt");

            if (width > 0 && hight > 0) {
                s2.SetImage(width, hight);
            }
            System.out.println("kmer=" + kmer);
            System.out.println("Initial string length =" + minlen);
            System.out.println("String length =" + seqlen);
            s2.Run();

            System.out.println("Shown repeated sequence is " + seqshow);
            if (flanksshow > 0) {
                System.out.println("Flanks around sequence is " + flanksshow);
            }

            long duration = (System.nanoTime() - startTime) / 1000000000;
            System.out.println("Time taken: " + duration + " seconds");
        } catch (IOException e) {
            System.out.println("Incorrect file name.");
        }
    }
}
