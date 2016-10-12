/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package reardist.GOCblastp.grimmgoctest;
import java.io.File;
import svasta.FileRW;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.text.DecimalFormat;

//todo: provjeri stdev i grimm values


/**
 *
 * Goal of the script:
 * to independently (from the previous code) calculate goc200 and produce files 
 * for grimm200 analysis.
 * The analysis will only take into account the data from the main chromosome, 
 * in cases where multiple genome elements are present.
 * @author Jelena
 */
public class GrimmGoc250 {
    

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        /*String ortdata = "C:/Users/Jelena/Desktop/transfer/newDataset/orthology/bidirectional_align_3/pblasting+dbs/bidirectionals_per_organism_bestHits&PercIdentityFilteredBranchOut/BranchOut/Enterococcus_faecium_DO_uid30627_longname_ortoLimit40";
        String positionData = "C:/Users/Jelena/Desktop/transfer/genbank/genomes/Bacteria/";
        String chromosomeList = "C:/Users/Jelena/Desktop/transfer/ListOfOrganismsAndGenElementsAll.txt";
        String grimfolder="C:/Users/Jelena/Desktop/transfer/newDataset/rearrangement_distances/grimm250New/Enterococcus_faecium/";
        */
        String ortdata = "C:\\wOrg\\writing_papers\\paper1\\newDataset\\orthologyTest\\bidirectionals_per_organism_bestHits&PercIdentityFilteredBranchOut\\";
        String positionData = "C:\\work\\genbank_05-10-2014\\ftp.ncbi.nih.gov\\genbank\\genomes\\Bacteria\\";
        String chromosomeList = "C:\\work\\genbank_05-10-2014\\ListOfOrganismsAndGenElementsAllTopology.txt";
        
        //for printing stuff out:
        String grimfoldy="C:\\wOrg\\writing_papers\\paper1\\newDataset\\rearrangement_distances_numOrthologs\\grimm100\\";
        String gocfolder="C:\\wOrg\\writing_papers\\paper1\\newDataset\\rearrangement_distances_numOrthologs\\goc100\\";
        new File(grimfoldy).mkdirs();
        new File(gocfolder).mkdirs();
        
        int numRepeats = 100;
        int numOrthologs = 100;
        Random random=new Random(3);
        File[] orthologdata=new File(ortdata).listFiles();
        
        
        for(File orts:orthologdata){
            if(orts.isDirectory() && !orts.getName().contains("Origins")){
                int num=0;
                System.out.println(orts.getName());
                String grimfolder=grimfoldy+orts.getName()+"/";
                new File(grimfolder).mkdirs();
                FileRW gocfile=new FileRW(gocfolder+orts.getName()+".txt","W");
                File[] pairs=orts.listFiles();
                for(File pair:pairs){
                    num++;
                    ArrayList<Ortholog> isOrdered=GrimmGoc250.sortOrthologs(pair, positionData, chromosomeList);
                    //isOrdered=GrimmgocFromBeginning3.readjustNameNumbers(isOrdered);
                    String grimsubfolder=grimfolder+orts.getName()+"\\"+orts.getName()+"_combination_"+num;
                    new File(grimsubfolder).mkdirs();
                    for(int i=0;i<100;i++){                        
                        ArrayList<Ortholog> orderedSubsample = GrimmGoc250.subsample(random, isOrdered, numOrthologs);
                        GrimmGoc250.toGrimmFile(orderedSubsample, new File(grimsubfolder+"\\" + pair.getName()+"_"+i));
                    }
                    //double goc=GrimmgocFromBeginning2.calculateGOC(isOrdered);
                    double[] gocs=GrimmGoc250.calculateAverageGOC(isOrdered, numRepeats, numOrthologs);
                    double averageGoc=GrimmGoc250.average(gocs);
                    double stdev=GrimmGoc250.stdev(gocs);
                    DecimalFormat df = new DecimalFormat("#.####");
                    gocfile.print(pair.getName()+"\t" + isOrdered.size()+"\t"+df.format(averageGoc)+"\t"+df.format(stdev)+"\n",false);
                }
                gocfile.print("", true);
            }
        }
    }
    
    
    public static double[] calculateAverageGOC(ArrayList<Ortholog> ordered, int numRepeats, int numOrthologs){
        
        double[] gocs=new double[numRepeats];
        Random random=new Random(3);
        for(int i=0;i<numRepeats;i++){
            ArrayList<Ortholog> orderedSubsample=GrimmGoc250.subsample(random, ordered, numOrthologs);
            double goc=GrimmGoc250.calculateGOC(orderedSubsample);
            gocs[i]=goc;
        }
        return gocs;
    }
    
    public static double average(double[] p) {
        double sum = 0;
        for (int i = 0; i < p.length; i++) {
            sum = sum + p[i];
        }
        double average = sum / p.length;
        return average;
    }
    
    public static double stdev(double[] p) {
        double sum = 0;  // sum of all the elements
        for (int i = 0; i < p.length; i++) {
            sum += p[i];
        }
        double mean=sum/p.length;
        double[] var=new double[p.length];
        for(int i = 0; i < p.length; i++){
            var[i] =Math.pow(mean-p[i], 2);
        }
        double varsum=0;
        for(int i = 0; i < p.length; i++){
            varsum+=var[i];
        }
        double stdev=Math.sqrt(varsum/(p.length-1));
        return stdev;
    }
    
    public static double calculateGOC(ArrayList<Ortholog> ordered){
        int[] pidPositions1=new int[ordered.size()];
        for(int i=0; i<ordered.size(); i++){
            pidPositions1[i]=ordered.get(i).getPid1position();
        }
        
        //in case the pidPositions for two genes map to the same spot, the output file will possibly contain a zero. Solution: exclude such genes -
        //the resulting genes-set will be smaller but only for a few genes
        pidPositions1=GrimmGoc250.extractUniqueInt(pidPositions1);
        ArrayList<Ortholog> or = new ArrayList<Ortholog>();
        for (int i : pidPositions1) {
            for (Ortholog o : ordered) {
                if (o.getPid1position() == i) {
                    or.add(o);
                    break;
                }
            }
        }
        ordered=or;
        
        
        //same thing for pidPositions2
        int[] pidPositions2=new int[ordered.size()];
        for(int i=0; i<ordered.size(); i++){
            pidPositions2[i]=ordered.get(i).getPid2position();
        }
        
        //in case the pidPositions for two genes map to the same spot, the output file will possibly contain a zero. Solution: exclude such genes -
        //the resulting genes-set will be smaller but only for a few genes
        pidPositions2=GrimmGoc250.extractUniqueInt(pidPositions2);
        or = new ArrayList<Ortholog>();
        for (int i : pidPositions2) {
            for (Ortholog o : ordered) {
                if (o.getPid2position() == i) {
                    or.add(o);
                    break;
                }
            }
        }
        ordered=or;
        
        if(pidPositions2.length<pidPositions1.length){
            pidPositions1=new int[pidPositions2.length];
            for (int i = 0; i < ordered.size(); i++) {
                pidPositions1[i] = ordered.get(i).getPid1position();
            }
        }
        
        if(pidPositions1.length!=250){
            System.out.println("Da, za goc - 2 pida bila na istom mjestu- " +" new dataset of size: " +pidPositions1.length+" "+pidPositions2.length);
        }
        if(pidPositions2.length!=250){
            System.out.println("Da, za goc - 2 pida na istom mjestu- " +" new dataset of size: " +pidPositions1.length+" "+pidPositions2.length);
        }
        
        Arrays.sort(pidPositions1);
        Arrays.sort(pidPositions2);
        
        int contiguous=0;
        int noncontiguous=0;
        
        for(int i=0; i<pidPositions1.length-1; i++){ //it shouldn't be minus 1!!!!!!!!!!!!!!!!!1
            //get the pair of gene positions that are contiguous in organism1
            int g1o1;
            int g2o1;
            try{
                g1o1=pidPositions1[i]; 
                g2o1=pidPositions1[i+1];
            }
            catch(IndexOutOfBoundsException e){ //due to the fact that bacterial chromosome is circular (well, it most often is :)
                g1o1=pidPositions1[i]; 
                g2o1=pidPositions1[0];
            };
            
            //identify orthologs they belong to and read their positions in organism2 
            Ortholog one=new Ortholog(); 
            Ortholog two=new Ortholog();
            for(Ortholog o:ordered){
                if(o.getPid1position()==g1o1){
                    one=o;
                }
                if(o.getPid1position()==g2o1){
                    two=o;
                }
            }
            int g1o2=one.getPid2position(); 
            int g2o2=two.getPid2position();
            
            if(g1o2==0 || g2o2==0){
                System.out.println("Sth wrong with recognizing correct orthologs; pid2 position is zero");
            }
            
            int index1=0;
            int index2=0;
            
            //are the positions of two orthologs contiguous in organism2
            for(int index=0;index<pidPositions2.length;index++){
                if(pidPositions2[index]==g1o2){
                    index1=index;
                }
                if(pidPositions2[index]==g2o2){
                    index2=index;
                }
            }
            
            if (index1==0 || index2==0){
                System.out.println("Sth wrong with recognizing correct orthologs; index is zero");
            }
            
            if(Math.abs((index1-index2))==1 || Math.abs((index1-index2))==(pidPositions2.length-1)){ //the two genes are near each other but are their orientations consistent?
                if((one.isPlusStrand1() && two.isPlusStrand1()) || (!one.isPlusStrand1() && !two.isPlusStrand1())){
                    if((one.isPlusStrand2() && two.isPlusStrand2()) || (!one.isPlusStrand2() && !two.isPlusStrand2())){
                        contiguous++;
                    }
                    else {
                        noncontiguous++;
                    }
                }
                else if((!one.isPlusStrand1() && two.isPlusStrand1()) || (one.isPlusStrand1() && !two.isPlusStrand1())){
                    if((!one.isPlusStrand2() && two.isPlusStrand2()) || (one.isPlusStrand2() && !two.isPlusStrand2())){
                        contiguous++;
                    }
                    else {
                        noncontiguous++;
                    }
                }
            }
            else{
                noncontiguous++;
            }
            
            
            
        }
        
        
        double goc=(double)contiguous/ordered.size();
        return goc;
    }
    

    
    
    public static void toGrimmFile(ArrayList<Ortholog> ordered, File file){
        int[] pidPositions1=new int[ordered.size()];
        for(int i=0; i<ordered.size(); i++){
            pidPositions1[i]=ordered.get(i).getPid1position();
        }
        
        //in case the pidPositions for two genes map to the same spot, the output file will possibly contain a zero. Solution: exclude such genes -
        //the resulting genes-set will be smaller but only for a few genes
        pidPositions1=GrimmGoc250.extractUniqueInt(pidPositions1);
        ArrayList<Ortholog> or = new ArrayList<Ortholog>();
        for (int i : pidPositions1) {
            for (Ortholog o : ordered) {
                if (o.getPid1position() == i) {
                    or.add(o);
                    break;
                }
            }
        }
        ordered=or;
        
        
        //same thing for pidPositions2
        int[] pidPositions2=new int[ordered.size()];
        for(int i=0; i<ordered.size(); i++){
            pidPositions2[i]=ordered.get(i).getPid2position();
        }
        
        //in case the pidPositions for two genes map to the same spot, the output file will possibly contain a zero. Solution: exclude such genes -
        //the resulting genes-set will be smaller but only for a few genes
        pidPositions2=GrimmGoc250.extractUniqueInt(pidPositions2);
        or = new ArrayList<Ortholog>();
        for (int i : pidPositions2) {
            for (Ortholog o : ordered) {
                if (o.getPid2position() == i) {
                    or.add(o);
                    break;
                }
            }
        }
        ordered=or;
        
        if(pidPositions2.length<pidPositions1.length){
            pidPositions1=new int[pidPositions2.length];
            for (int i = 0; i < ordered.size(); i++) {
                pidPositions1[i] = ordered.get(i).getPid1position();
            }
        }
        
        
        if(pidPositions1.length!=250){
            System.out.println("Da, 2 pida bila na istom mjestu za file1: " +file.getAbsolutePath()+" new dataset of size: " +pidPositions1.length+" "+pidPositions2.length);
        }
        if(pidPositions2.length!=250){
            System.out.println("Da, 2 pida na istom mjestu za file2: " +file.getAbsolutePath()+" new dataset of size: " +pidPositions1.length+" "+pidPositions2.length);
        }
        
        
        
        Arrays.sort(pidPositions1);
        Arrays.sort(pidPositions2);
        
        int[] toprint1=new int[pidPositions1.length];
        int[] toprint2=new int[pidPositions2.length];
        
        for(int i=0; i<pidPositions1.length-1; i++){
            //get the pair of gene positions that are contiguous in organism1
            int g1o1=pidPositions1[i]; 
            int g2o1=pidPositions1[i+1];
            
            
            
            //identify orthologs they belong to and read their positions in organism2 
            Ortholog one=new Ortholog(); 
            Ortholog two=new Ortholog();
            for(Ortholog o:ordered){
                if(o.getPid1position()==g1o1){
                    one=o;
                }
                if(o.getPid1position()==g2o1){
                    two=o;
                }
            }
            int g1o2=one.getPid2position(); 
            int g2o2=two.getPid2position();

            //while you are at it, put the orthologs in the toprint1 array
            if(one.isPlusStrand1()){toprint1[i]=(i+1);}
            else{toprint1[i]=-(i+1);}
            if(two.isPlusStrand1()){toprint1[i+1]=(i+1)+1;}
            else{toprint1[i+1]=-(i+1+1);}
            
            
            //are the positions of two orthologs contiguous in organism2
            int index1=0;
            int index2=0;
            for(int index=0;index<pidPositions2.length;index++){
                if(pidPositions2[index]==g1o2){
                    index1=index;
                }
                if(pidPositions2[index]==g2o2){
                    index2=index;
                }
            }
            
            //put the orthologs in the toprint2 array
            if(one.isPlusStrand2()){toprint2[index1]=i+1;}
            else{toprint2[index1]=-(i+1);}
            if(two.isPlusStrand2()){toprint2[index2]=i+1+1;}
            else{toprint2[index2]=-(i+1+1);}
        }
        
        //check: what's up with the comparison of the last and first ortholog in the array?
        //check: is it relevant which are minus and which plus genes in the second genome - i.e. do things change if i reverse it from the current calls

        for(int i=0;i<toprint1.length;i++){
            if(toprint1[i]==0){
                System.out.println("zero present in file "+file.getAbsolutePath());
            }
        }
        for(int i=0;i<toprint2.length;i++){
            if(toprint2[i]==0){
                System.out.println("zero present in file "+file.getAbsolutePath());
            }
        }
        
        
        FileRW f=new FileRW(file.getAbsolutePath(), "W");
        f.print(">"+file.getName().split(";")[0]+"\n", false);
        for(int o:toprint1){
            f.print(o+" ", false);
        }
        f.print("\n>"+file.getName().split(";")[1]+"\n", false);
        for(int o:toprint2){
            f.print(o+" ", false);
        }
        f.print("", true);
    }
    
    
    public static ArrayList<Ortholog> subsample(Random r, ArrayList<Ortholog> a, int numOrthologs){
        ArrayList<Ortholog> subsample=new ArrayList<Ortholog>();
        
        if (a.size() <= (numOrthologs)) {
            subsample = a;
            return subsample;                    
        }
        
        int[] indices=new int[numOrthologs]; 
        ArrayList<Integer> potentialIndices=new ArrayList<Integer>(); //contains actual indices of the elements in the arraylist a, the indices will be subsampled
        
        for(int i=0; i<a.size(); i++){
            potentialIndices.add(i);
        }
        
        for(int i=0; i<numOrthologs; i++){
            int k=r.nextInt(potentialIndices.size()); //check the method, does it include zero; yes it does
            indices[i]=potentialIndices.get(k); 
            //indices[i]=k;
            potentialIndices.remove(k); //sampling without replacement
        }
        
        Arrays.sort(indices); //the purpose of this is for the sublist elements to be in the same order as in the incoming list
       
        for(int i=0; i<indices.length; i++){
            subsample.add(a.get(indices[i]));
        }
        
        return subsample;
    }
    
    public static ArrayList<Ortholog> sortOrthologs(File pair, String positionData, String chromosomeList){ //ArrayList<Ortholog>
        String[] orgs=pair.getName().split(";");
        File org1=new File(positionData+orgs[0]);
        String chr1="";
        File org2=new File(positionData+orgs[1]);
        String chr2="";
        
        //find chromosome names for the two organisms
        FileRW f=new FileRW(chromosomeList, "R");
        String line=f.readLine();
        while((line!=null)){
            if(line.contains(orgs[0])){
                int maxlength=0;
                line=f.readLine();
                while (!line.contains(">") && line.length()>0 && line!=null){
                    int chrlength=Integer.parseInt(line.split("\\t")[2]);
                    if(chrlength>maxlength){
                        chr1=line.split("\\t")[0];
                        maxlength=chrlength;
                    }
                    line=f.readLine();
                } 
            }
            line=f.readLine();
        }
        
        f=new FileRW(chromosomeList, "R");
        line=f.readLine();
        while((line!=null)){
            if(line.contains(orgs[1]) && line!=null){
                int maxlength=0;
                line=f.readLine();
                try{
                while (line!=null && !line.contains(">") && line.length()>0){
                    int chrlength=Integer.parseInt(line.split("\\t")[2]);
                    if(chrlength>maxlength){
                        chr2=line.split("\\t")[0];
                        maxlength=chrlength;
                    }
                    line=f.readLine();
                }
                }catch(Exception e){
                    System.out.println();
                }
            }
            line=f.readLine();
        }
        
        String c1= org1.getAbsolutePath()+"\\"+chr1+".ptt";
        String c2= org2.getAbsolutePath()+"\\"+chr2+".ptt";

        
        ArrayList<Ortholog> o=new ArrayList<Ortholog>();
        //get the PID data from the ortholog files into the arraylist
        f=new FileRW(pair.getAbsolutePath(), "R");
        
        //read in the data on the BBHs
        while((line=f.readLine())!=null){
            String[] data=line.split("\\|");
            Ortholog ort=new Ortholog();
            ort.setPid1(data[1]);
            ort.setPid2(data[5]);
            o.add(ort);
    
        }
        
        f=new FileRW(c1, "R"); //get the positional data on orthologs from the genbank .ptt files
        while((line=f.readLine())!=null){
            for(Ortholog ort:o){
                if(line.contains(ort.getPid1())){
                    ort.setPtt1(line);
                    String[] s=line.split("\\s++");
                    if(s[1].contains("+")){
                        ort.setPlusStrand1(true);
                    }
                    String[] k=s[0].split("\\.\\.");
                    if(ort.isPlusStrand1()){
                        ort.setPid1position(Integer.parseInt(k[0]));
                    }
                    else{
                        ort.setPid1position(Integer.parseInt(k[1]));
                    }
                    
                }
            }
        }
        
        f=new FileRW(c2, "R"); //get the positional data on orthologs from the genbank .ptt files
        while((line=f.readLine())!=null){
            for(Ortholog ort:o){
                if(line.contains(ort.getPid2())){
                    ort.setPtt2(line);
                    String[] s=line.split("\\s++");
                    if(s[1].contains("+")){
                        ort.setPlusStrand2(true);
                    }
                    String[] k=s[0].split("\\.\\.");
                    if(ort.isPlusStrand2()){
                        ort.setPid2position(Integer.parseInt(k[0]));
                    }
                    else{
                        ort.setPid2position(Integer.parseInt(k[1]));   //check results
                    }
                }
            }
        }
        
        ArrayList<Ortholog> toreturn=new ArrayList<Ortholog>();
        for(Ortholog ort:o){
            if(ort.getPid2position()!=0 && ort.getPid1position()!=0){ //ortPosition equals 0 in cases where the gene wasn't found on the given genome element (or very unlikely, if the gene actually starts at position 0)
                toreturn.add(ort);
            }
            else{
                //System.out.println();
            }
        }
        return toreturn;
    }   
    
    
    public static int[] extractUniqueInt(int[] inty){
        ArrayList<Integer> intt=new ArrayList<Integer>();
        for(int i:inty){
            if(!intt.contains(i)){
                intt.add(i);
            }
        }
        
        int[] toreturn=new int[intt.size()];
        for(int i=0;i<intt.size();i++){
            toreturn[i]=intt.get(i);
        }
        return toreturn;
        
    }
    
}
