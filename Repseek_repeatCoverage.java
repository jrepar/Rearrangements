/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package repeats.repseek;

import java.io.File;
import java.util.ArrayList;
import svasta.FileRW;

/**
 * * Goal of script:
 * Graph that shows which part of chromosomes are most likely to be affected by genome rearrangements
 * Rearrangement is affecting the region between two repeats, the shorter one (when taking account the circular chromosome)
 * 
 * 
 * @author jrepar
 */
public class Repseek_repeatCoverage {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        File out=new File("C:\\J\\Origins\\repeatitiveness\\repseek_analysisKnownOris_repeatRegionsCoverageInversionsWeighed_seedOnly_p0.001_chromosomes\\");
        //File out=new File("C:\\J\\Origins\\repeatitiveness\\repseek_analysis_repeatRegionsCoverageDirect_seedOnly_p0.001_chromosomes\\");
        out=new File("E:\\repseekAllAnalysis\\repseek_analysis_repeatRegionsCoverageInversions_seedOnly_p0.001_chromosome\\");
        
        File foldy=new File("E:\\repseek_results_treshold10_p0.001_chromosomesOrigins\\");
        foldy=new File("E:\\repseekAll\\repseek_results_p0.001_chromosomeLargestGelement\\");
        
        out.mkdirs();
        String genomes="E:\\gen_refseq_NCBI_genomes\\ListOfOrganismsAndGenElementsAllTopology_DraftMarked.txt";
        //String genomes="E:\\genbank_05-10-2014\\ListOfOrganismsAndGenElementsAllTopology.txt";
       
        boolean weighed=false;
        String invOrDirect="i";
        
        
        //read in the refseq/genbank data to filter the genome elements
        ArrayList<String> ao = new ArrayList<String>();
        ArrayList<ArrayList<String>> bo = new ArrayList<ArrayList<String>>();
        ArrayList<String> joj = new ArrayList<String>();
        FileRW refList = new FileRW(genomes, "R");
        String line = refList.readLine();
        ao.add(line);
        while ((line = refList.readLine()) != null) {
            if (line.startsWith(">")) {
                ao.add(line);
                bo.add(joj);
                joj = new ArrayList<String>();
            } else {
                joj.add(line);
            }
        }
        bo.add(joj);
        
        if(foldy.isDirectory()){
            File[] orgs = foldy.listFiles();
            for (File focus : orgs) {
                File out2=new File(out.getAbsolutePath()+"\\"+focus.getName());
                if(focus.getName().contains("doriC") ){
                    out2.mkdirs();
                    File[] org=focus.listFiles();
                    for(File c:org){
                        if (c.getName().startsWith("out")) {
                            //find the genome length of the largest gelement
                            String o=c.getName().substring(3).split("\\.")[0];
                            ArrayList<String> gelements=new ArrayList<String>();
                            for(int i=0;i<ao.size();i++){
                                String oo=ao.get(i);
                                if(oo.contains(o)){
                                    gelements=bo.get(i);
                                }
                            }
                            int max=0;
                            int indexMax=342;
                            for(int k=0;k<gelements.size();k++){
                                String gelement=gelements.get(k);
                                int maxy=Integer.parseInt(gelement.split("\t")[2]);
                                if(maxy>max){
                                    max=maxy;
                                    indexMax=k;
                                }
                            }

                            int[] genome=new int[max];

                            FileRW r = new FileRW(c.getAbsolutePath(), "R");
                            while ((line = r.readLine()) != null) {
                                if (line.startsWith(invOrDirect)) { //direct repeats or inversions!!!!!!!!!!!!!!
                                    int one=Integer.parseInt(line.split("\t")[1]);
                                    int two=Integer.parseInt(line.split("\t")[2]);

                                    int miny=Math.min(two, one);
                                    int maxy=Math.max(two, one);

                                    int gsize1=maxy-miny;
                                    int gsize2=max-gsize1;
                                    if(!weighed){
                                        if(gsize1>gsize2){
                                            for(int u=maxy;u<genome.length;u++){
                                                genome[u]++;
                                            }
                                            for(int u=0;u<miny;u++){
                                                genome[u]++;
                                            }
                                        }
                                        else{
                                            for(int u=miny; u<maxy;u++){
                                                genome[u]++;
                                            }
                                        }
                                    }
                                    else{
                                        int weight=Integer.parseInt(line.split("\t")[3]);
                                        if (gsize1 > gsize2) {
                                            for (int u = maxy; u < genome.length; u++) {
                                                genome[u]=genome[u]+weight;
                                            }
                                            for (int u = 0; u < miny; u++) {
                                                genome[u]=genome[u]+weight;
                                            }
                                        } else {
                                            for (int u = miny; u < maxy; u++) {
                                                genome[u]=genome[u]+weight;
                                            }
                                        }
                                    }
                                }
                            }
                            FileRW w = new FileRW(out2.getAbsolutePath() +"\\"+c.getName().substring(3).split("\\.")[0], "W");
                            w.print("GenomeLocation\tRearrCoverage\n", false);
                            for(int i=0;i<genome.length;i++){
                                w.print(i+"\t"+genome[i]+"\n", false);
                            }
                            w.print("",true);
                        }
                    }
                }
                
            }    
        }
    }
    
}
