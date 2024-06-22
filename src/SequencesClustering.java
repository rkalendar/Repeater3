import java.util.HashMap;

public final class SequencesClustering {

    public SequencesClustering(String[] primers, int sim, boolean t, int kmer, int minlenseq) {
        this.seq = primers;
        if (sim < 60) {
            sim = 60;
        }
        if (sim > 90) {
            sim = 90;
        }
        if (kmer < 12) {
            kmer = 12;
        }
        if (kmer > 99) {
            kmer = 99;
        }
        this.sim = sim;
        nseq = primers.length;
        ncl = Clustering(kmer, minlenseq);
    }

    public int[] Result() {
        return cx;
    }

    private int Clustering(int kmer, int minlenseq) {
        int n = 0;
        String s;
        String a;
        cx = new int[nseq];      // cluster
        int[] lx = new int[nseq];

        HashMap<String, int[]> map = new HashMap<>();

        for (int j = 0; j < nseq; j++) {
            lx[j] = seq[j].length();
            if (lx[j] > minlenseq) {
                a = dna.ComplementDNA(seq[j]);
                for (int i = kmer - 1; i < lx[j]; i++) {
                    s = seq[j].substring(i - kmer + 1, i + 1);
                    if (map.containsKey(s)) {
                        int[] t = map.get(s);
                        t[j]++;
                        map.put(s, t);
                    } else {
                        int[] t = new int[nseq];
                        t[j] = 1;
                        map.put(s, t);
                    }
                    s = a.substring(i - kmer + 1, i + 1);
                    if (map.containsKey(s)) {
                        int[] t = map.get(s);
                        t[j]++;
                        map.put(s, t);
                    } else {
                        int[] t = new int[nseq];
                        t[j] = 1;
                        map.put(s, t);
                    }
                }
            }
        }

        //clustering
        for (int j = 0; j < nseq; j++) {
            if (lx[j] > minlenseq) {
                int[] x = new int[nseq];
                a = dna.ComplementDNA(seq[j]);
                for (int i = kmer - 1; i < lx[j]; i++) {
                    s = seq[j].substring(i - kmer + 1, i + 1);
                    int[] t = map.get(s);
                    for (int h = 0; h < nseq; h++) {
                        if (t[h] > 0) {
                            x[h]++;
                        }
                    }
                    s = a.substring(i - kmer + 1, i + 1);
                    t = map.get(s);
                    for (int h = 0; h < nseq; h++) {
                        if (t[h] > 0) {
                            x[h]++;
                        }
                    }
                }

                for (int h = 0; h < nseq; h++) {
                    if (x[h] > 0) {
                        int m = (100 * x[h]) / x[j];
                        if (m > sim) {
                            if (cx[j] == 0) {
                                if (cx[h] > 0) {
                                    cx[j] = cx[h];
                                } else {
                                    n++;
                                    cx[j] = n;
                                    cx[h] = n;
                                }
                            } else {
                                if (cx[h] > 0) {
                                    if (cx[h] < cx[j]) {
                                        cx[j] = cx[h];
                                    } else {
                                        cx[h] = cx[j];
                                    }
                                } else {
                                    cx[h] = cx[j];
                                }

                            }

                        }
                    }
                }
            }
        }
        return n;
    }

    public int getNcl() {
        return ncl;
    }
    private final int nseq;
    private final int sim;
    private int ncl;
    private int[] cx;
    private final String seq[];
}
