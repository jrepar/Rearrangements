/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package operons;
import java.io.File;
import svasta.FileRW;
import java.util.ArrayList;
import java.text.DecimalFormat;

/**
 *
 * @author jrepar
 */
public class OperonsQuantifyTrainSize {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        String genbank="C:\\J\\genbank_05-10-2014\\ftp.ncbi.nih.gov\\genbank\\genomes\\Bacteria\\";
        String orgList="C:\\J\\paper1\\newDataset\\phylogeny\\16srrna\\BranchoutApproach\\longnames_20\\TheChosenOnesFinal\\";
        String genomeList="C:\\J\\genbank_05-10-2014\\ListOfOrganismsAndGenElementsAllTopology.txt";
        String out="C:\\J\\paper1\\newDataset\\operons\\trainSize\\";
        
        
        FileRW refList = new FileRW(genomeList, "R");
        ArrayList<String> ab = new ArrayList<String>(); // will contain organism names
        ArrayList<ArrayList<String>> bo = new ArrayList<ArrayList<String>>(); //will contain data on genome elements corresponding to organism names in ArrayList ab
        ArrayList<String> joj = new ArrayList<String>();
        String line = refList.readLine();
        ab.add(line);
        while ((line = refList.readLine()) != null) {
            if (line.startsWith(">")) {
                ab.add(line);
                bo.add(joj);
                joj = new ArrayList<String>();
            } else {
                joj.add(line);
            }
        }
        bo.add(joj);
        
        
        File[] orgs=new File(orgList).listFiles();
        for(File o:orgs){
                FileRW w=new FileRW(out+o.getName(),"W");
                w.print("ORGANISM\tSizeOfGenome\tNumOfGenes\tAverageTrainLengthAll\tStdevTrainLengthAll\tNumOfTrainsAll\tAverageTrainLengthOnes\tStdevTrainLengthOnes\tNumOfTrainsOnes\tAverageTrainLengthOperons\tStdevTrainLengthOperons\tNumOfTrainsOperons\n", false);
                FileRW r=new FileRW(o.getAbsolutePath(),"R");
                while((line=r.readLine())!=null){
                    String organism=line.trim();
                    //find the name of the largest genome element of the organism
                    for(int i=0;i<ab.size();i++){
                        String str=ab.get(i);
                        if(str.endsWith(organism)){
                            
                            ArrayList<String> gelements=bo.get(i);
                            String maxG="jač";
                            int max=0;
                            for(int k=0;k<gelements.size();k++){
                                String g=gelements.get(k);
                                int gsize=Integer.parseInt(g.split("\t")[2]);
                                if(gsize>max){
                                    maxG=g.split("\t")[0];
                                    max=gsize;
                                }
                                
                            }
                            File bact=new File(genbank+organism);
                            File[] files=bact.listFiles();
                            ArrayList<Integer> ops=new ArrayList<Integer>();
                            int numGenes=0;
                            int size=0;
                            String direction="sfkča";
                            for(File f:files){
                                if(f.getName().endsWith("ptt") && f.getName().contains(maxG)){
                                    FileRW rd=new FileRW(f.getAbsolutePath(),"R");
                                    String gene=rd.readLine();
                                    gene=rd.readLine();
                                    gene=rd.readLine();
                                    while((gene=rd.readLine())!=null){
                                        if(line.length()>0){
                                            numGenes++;
                                            if(gene.split("\t")[1].equalsIgnoreCase(direction)){
                                                size++;
                                            }
                                            else if(!gene.split("\t")[1].equalsIgnoreCase(direction) && !direction.equalsIgnoreCase("sfkča")){
                                                ops.add(size);
                                                direction=gene.split("\t")[1];
                                                size=1;
                                            }
                                            else if(direction.equalsIgnoreCase("sfkča")){
                                                direction=gene.split("\t")[1];
                                                size=1;
                                            }
                                        }
                                    }
                                    ops.add(size);
                                }
                            }
                            ArrayList<Integer> ones=new ArrayList<Integer>();
                            ArrayList<Integer> moreThanOnes=new ArrayList<Integer>();
                            for(Integer inty:ops){
                                if(inty==1){
                                    ones.add(inty);
                                }
                                else{
                                    moreThanOnes.add(inty);
                                }
                            }
                            
                            double average=OperonsQuantifyTrainSize.average(ops);
                            double stdev=OperonsQuantifyTrainSize.stdev(ops);
                            DecimalFormat df = new DecimalFormat("#.##");
                            //w.print("ORGANISM\tSizeOfGenome\tNumOfGenes\tAverageTrainLengthAll\tStdevTrainLengthAll\tNumOfTrainsAll\tAverageTrainLengthOnes\tStdevTrainLengthOnes\tNumOfTrainsOnes\tAverageTrainLengthOperons\tStdevTrainLengthOperons\tNumOfTrainsOperons\n", false);
                
                            w.print(organism+"\t"+Integer.toString(max)+"\t"+numGenes+"\t"+df.format(average)+"\t"+df.format(stdev)+"\t"+ops.size()+"\t"+df.format(OperonsQuantifyTrainSize.average(ones))+"\t"+df.format(OperonsQuantifyTrainSize.stdev(ones))+"\t"+ones.size()+"\t"+df.format(OperonsQuantifyTrainSize.average(moreThanOnes))+"\t"+df.format(OperonsQuantifyTrainSize.stdev(moreThanOnes))+"\t"+moreThanOnes.size()+"\n", false);
                        }
                    }
                }
                w.print("", true);
            
        }
    }
    
    public static double average(ArrayList<Integer> p) {
        double sum = 0;
        for (int i = 0; i < p.size(); i++) {
            sum = sum + p.get(i);
        }
        double average = (double) sum / p.size();
        return average;
    }
    
    public static double stdev(ArrayList<Integer> p) {
        double sum = 0;  // sum of all the elements
        for (int i = 0; i < p.size(); i++) {
            sum += p.get(i);
        }
        double mean=sum/p.size();
        double[] var=new double[p.size()];
        for(int i = 0; i < p.size(); i++){
            var[i] =Math.pow(mean-p.get(i), 2);
        }
        double varsum=0;
        for(int i = 0; i < p.size(); i++){
            varsum+=var[i];
        }
        double stdev=Math.sqrt((double)varsum/(p.size()-1));
        return stdev;
    }
}
