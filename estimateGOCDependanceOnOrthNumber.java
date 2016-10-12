/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package reardist.GOCblastp;
import svasta.FileRW;
import java.util.ArrayList;
import java.io.File;
import java.text.DecimalFormat;
import java.util.Random;
import java.util.Arrays;

/**
 * Goal:
 * Choose a set of organisms and on them do a simulation:
 * choose a starting number of orthologs and a number you are going to increase it with (eg 20,20).
 * For each organism pair, choose 20 (or the chosen start num) orthologs randomly,
 * and calculate GOC based on these orthologs.
 * Also, choose a number of repeats of this experiment that happens at a certain number of orthologs (eg 1000).
 * Report average and stdev at each point (i.e. at each number of orthologs).
 *
 *
 * @author Jelena
 */
public class estimateGOCDependanceOnOrthNumber {


    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        int startnum2=100;
        int increasenum=100;
        int repeatExperiment=100;
        Random random=new Random();
        System.out.println("### NumOfResamplingsPerNumOfOrthologs: "+repeatExperiment);
        System.out.println(">organism");
        System.out.println("NumOfOrthologs\tAverageGOC\tGOCstdev");
        //System.out.println("contiguousGenes noncontiguousGenes  numOfOrthologs  GOC");
        String path="C:\\Users\\jrepar\\Desktop\\work\\Orthology\\bidirectional_align_3\\pblasting+dbs\\bidirectionals_per_organism_bestHits&PercIdentityFilteredBranchOut\\T_radiovictrix\\";
        path="C:\\J\\paper1\\newDataset\\orthologyFinal\\bidirectionals_per_organism_bestHits&PercIdentityFilteredBranchOut\\T_gammatolerans\\";
        File file=new File(path);
        String path2="C:\\Users\\jrepar\\Desktop\\work\\genbank_05-10-2014\\ftp.ncbi.nih.gov\\genbank\\genomes\\Bacteria\\";
        path2="C:\\J\\genbank_05-10-2014\\ftp.ncbi.nih.gov\\genbank\\genomes\\Bacteria\\";
        if(file.isDirectory()){ //file should contain ortholog data of organism pairs
            File[] fileList = file.listFiles();
            for(File f:fileList){
                if(!f.getName().contains("Anaeromyxobacter") && !f.getName().contains("Syntrophobacter") && !f.getName().contains("Desulfobacca")){ //&& f.getName().contains("Truepera") && f.getName().contains("maricope")
                    System.out.println(">"+f.getName()+"\t");
                    //System.out.println();
                    ArrayList<String> orthologs=new ArrayList<String>();
                    if(!f.isDirectory()){
                        //read in pairwise ortholog data
                        FileRW it = new FileRW(f.getAbsolutePath(), "R");
                        String line = "";
                        while ((line = it.readLine()) != null) {
                            if(line.length()>0){
                                orthologs.add(line);
                            }
                        }
                        String[] orgs=f.getName().split(";");
                        File pathorg1=new File(path2+orgs[0]);
                        File pathorg2=new File(path2+orgs[1]);
                        
                        if (orgs.length == 3) {
                            if (pathorg1.exists()) {
                                pathorg2 = new File(path2+ orgs[1] + "-" + orgs[2]);
                            } else {
                                pathorg1 = new File(path2 + orgs[0] + "-" + orgs[1]);
                                pathorg2 = new File(path2 + orgs[2]);
                            }
                        }
                        else if (f.getName().contains("Shewanella_W3-18-1_uid13902")) {
                            String[] g = f.getName().split("Shewanella_W3-18-1_uid13902");
                            if (g.length > 1 && g[1].startsWith("-")) {
                                pathorg1 = new File(path2 + "Shewanella_W3-18-1_uid13902");
                                pathorg2 = new File(path2 + g[1].substring(1));
                            } else {
                                pathorg1 = new File(path2 + g[0].substring(0, g[0].length() - 1));
                                pathorg2 = new File(path2 + "Shewanella_W3-18-1_uid13902");

                            }
                        }
                        else if (orgs.length == 4) {
                            pathorg1 = new File(path2 + orgs[0] + "-" + orgs[1]);
                            pathorg2 = new File(path2 + orgs[2] + "-" + orgs[3]);
                        } else if (!(orgs.length == 2)) {
                            System.out.println("Warning!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
                        }

                        ArrayList<ArrayList<String>> org1=new ArrayList<ArrayList<String>>();
                        ArrayList<ArrayList<String>> org2=new ArrayList<ArrayList<String>>();
                        //read in gene position data from ncbi files for organism1
                        if(pathorg1.isDirectory()){
                            File[] filei = pathorg1.listFiles();
                            for(File fil:filei){
                                //ArrayList<String> a=new ArrayList<String>();
                                if(fil.getName().endsWith("ptt")){
                                    ArrayList<String> a=new ArrayList<String>();
                                    FileRW t=new FileRW(fil.getAbsolutePath(), "R");
                                    for(int o=0;o<3;o++){
                                        line=t.readLine();
                                    }
                                    line="";
                                    while ((line = t.readLine()) != null) {
                                        if(line.length()>0 && line.contains("..")){
                                            //use only genes you have identified orthologs for
                                            String[] k=line.split("\t");
                                            for(String ort:orthologs){
                                                if(ort.contains(k[3])){
                                                    a.add(line);
                                                }
                                            }
                                        }
                                    }
                                    org1.add(a);
                                }
                                //org1.add(a);
                            }
                        }
                        //read in gene position data from ncbi files for organism2
                        if(pathorg2.isDirectory()){
                            File[] filei = pathorg2.listFiles();
                            for(File fil:filei){
                                //ArrayList<String> a=new ArrayList<String>();
                                if(fil.getName().endsWith("ptt")){
                                    ArrayList<String> a=new ArrayList<String>();
                                    FileRW t=new FileRW(fil.getAbsolutePath(), "R");
                                    for(int o=0;o<3;o++){
                                        line=t.readLine();
                                    }
                                    line="";
                                    while ((line = t.readLine()) != null) {
                                        if(line.length()>0 && line.contains("..")){
                                            //use only genes you have identified orthologs for
                                            String[] k=line.split("\t");
                                            for(String ort:orthologs){
                                                if(ort.contains(k[3])){
                                                    a.add(line);
                                                }
                                            }
                                        }
                                    }
                                    org2.add(a);
                                }
                            }
                        }



                        //Identify biggest geneleemnts of org1 and org2, i.e. their
                        //indices in the arraylist org1 and org2
                        int ge1=0;
                        int ge2=0;
                        int ge1size=0;
                        int ge2size=0;
                        for(int h=0;h<org1.size();h++){
                            if(org1.get(h).size()>ge1size){
                                ge1=h;
                                ge1size=org1.get(h).size();
                            }
                        }
                        for(int h=0;h<org2.size();h++){
                            if(org2.get(h).size()>ge2size){
                                ge2=h;
                                ge2size=org2.get(h).size();
                            }
                        }

                        //Each item in the ArrayList orthologs is data for one ortholog
                        //Choose randomly startnum number of orthologs and calculate GOC for them.
                        //repeat this repeatExperiment times
                        //Increase the startnum by increasenum and repeat.
                        //which organisms am I calculating GOC for
                        int startnum=startnum2;
                        int smaller;
                        if(org1.get(ge1).size()<org2.get(ge2).size()){
                            smaller=org1.get(ge1).size();
                        }
                        else{
                            smaller=org2.get(ge2).size();
                        }
                        while(startnum<=(smaller-500)){
                            //System.out.println("Num of orthologs: " + startnum);
                            System.out.print(startnum+"\t");
                            double[] collectGOC=new double[repeatExperiment];
                            for(int u=0; u<repeatExperiment;u++){
                                int contig=0;
                                int noncontig=0;
                                //I only wish to take the GOC values for one genome
                                //element. Because, the purpose of this exercise is
                                //to see if I will be able to accurately measure
                                //rearrangements on separate genome elements, including
                                //the small ones, with GOC. For this, the longest genome element
                                //will be the most useful one if I wish to see the
                                //influence of the number of orthologs on calculated GOC.
                                //Therefore, in choosing the orthologs, I must be
                                //careful to only take those belonging to the chromosome
                                //(or bigger chromosome if there are two)
                                String[] reducedOrthologs=new String[startnum];
                                String[] toBeOrderedOrthologPTTsOrg1=new String[startnum];
                                String[] toBeOrderedOrthologPTTsOrg2=new String[startnum];
                                ArrayList<String> workingOrts=(ArrayList<String>)orthologs.clone();
                                for(int howmany=0; howmany<startnum; howmany++){
                                    int randy=random.nextInt(workingOrts.size());
                                    String ort=workingOrts.get(randy);
                                    workingOrts.remove(randy);
                                    /*ArrayList<String> wo=new ArrayList<String>();
                                    for(String h: workingOrts){
                                        if(workingOrts.indexOf(h)!=randy){
                                            wo.add(h);
                                        }
                                    }
                                    workingOrts=wo;*/
                                    String[] h = ort.split("\t")[1].split("\\|");
                                    //the two pids are now in h[1] and h[5] positions
                                    int totake=0;
                                    for(String or:org1.get(ge1)){
                                        if(or.contains(h[1]) || or.contains(h[5])){
                                            totake++;
                                            continue;
                                        }
                                    }
                                    for(String or:org2.get(ge2)){
                                        if(or.contains(h[1]) || or.contains(h[5])){
                                            totake++;
                                            continue;
                                        }
                                    }
                                    int hasit = 0; //for sampling without replacement
                                    for (String o : reducedOrthologs) {
                                        if (o != null && o.equalsIgnoreCase(ort)) {
                                            hasit++;
                                        }
                                    }
                                    if(totake==2 && hasit==0) {
                                        reducedOrthologs[howmany] = ort;

                                        for (String or : org1.get(ge1)) {
                                            if (or.contains(h[1]) || or.contains(h[5])) {
                                                toBeOrderedOrthologPTTsOrg1[howmany] = or;
                                                continue;
                                            }
                                        }
                                        for (String or : org2.get(ge2)) {
                                            if (or.contains(h[1]) || or.contains(h[5])) {
                                                toBeOrderedOrthologPTTsOrg2[howmany] = or;
                                                continue;
                                            }
                                        }
                                    }
                                    else{
                                        howmany--;
                                    }

                                }

                                //order orthologs by their position on the chromosome in both organism1 and 2    //!!!!!!!!!!!!!!!!!!!!!!!!!do i have to choose the orthologs that belong ONLY to the chromosome???
                                int[] orderedOrthologPTTsOrg1=new int[startnum];
                                int[] orderedOrthologPTTsOrg2=new int[startnum];
                                for(int r=0; r<toBeOrderedOrthologPTTsOrg1.length;r++){
                                    String[] uh=toBeOrderedOrthologPTTsOrg1[r].split("\t");
                                    String start = uh[0].split("\\.\\.")[0];
                                    orderedOrthologPTTsOrg1[r]=Integer.valueOf(start);
                                }
                                for(int r=0; r<toBeOrderedOrthologPTTsOrg2.length;r++){
                                    String[] uh=toBeOrderedOrthologPTTsOrg2[r].split("\t");
                                    String start = uh[0].split("\\.\\.")[0];
                                    orderedOrthologPTTsOrg2[r]=Integer.valueOf(start);
                                }

                                Arrays.sort(orderedOrthologPTTsOrg1);
                                Arrays.sort(orderedOrthologPTTsOrg2);

                                //Start comparisons-which genes are contiguous, and which are not
                                int startx;
                                int starty;
                                for (int i = 0; i < orderedOrthologPTTsOrg1.length; i++) {
                                    if (i == (orderedOrthologPTTsOrg1.length - 1)) { //because the chromosome is a circle
                                        startx= orderedOrthologPTTsOrg1[i];
                                        starty = orderedOrthologPTTsOrg1[0];
                                    } else {
                                        startx = orderedOrthologPTTsOrg1[i];
                                        starty = orderedOrthologPTTsOrg1[i+1];
                                    }
                                    //which genes are these in organism1?
                                    int genexindex=0;
                                    int geneyindex=0;
                                    for(int j=0;j<toBeOrderedOrthologPTTsOrg1.length; j++){
                                        String[] uh=toBeOrderedOrthologPTTsOrg1[j].split("\t");
                                        if(uh[0].contains(Integer.toString(startx))){
                                            genexindex=j;
                                        }
                                        if(uh[0].contains(Integer.toString(starty))){
                                            geneyindex=j;
                                        }
                                    }
                                    //which genes are these in organism2? they are
                                    //the genes at the positions genexindex and geneyindex
                                    //in the array toBeOrderedOrthologPTTsOrg2
                                    int startx2=Integer.valueOf(toBeOrderedOrthologPTTsOrg2[genexindex].split("\t")[0].split("\\.\\.")[0]);
                                    int starty2=Integer.valueOf(toBeOrderedOrthologPTTsOrg2[geneyindex].split("\t")[0].split("\\.\\.")[0]);
                                    //and what is the relative position of these genes in organism2
                                    int index=0;
                                    int indey=0;
                                    for(int ik = 0; ik < orderedOrthologPTTsOrg2.length; ik++){
                                        if(orderedOrthologPTTsOrg2[ik]==startx2){
                                            index=ik;
                                        }
                                        if(orderedOrthologPTTsOrg2[ik]==starty2){
                                            indey=ik;
                                        }
                                    }
                                    String orix1 = toBeOrderedOrthologPTTsOrg1[genexindex].split("\t")[1];
                                    String orix2 = toBeOrderedOrthologPTTsOrg2[genexindex].split("\t")[1];
                                    String oriy1 = toBeOrderedOrthologPTTsOrg1[geneyindex].split("\t")[1];
                                    String oriy2 = toBeOrderedOrthologPTTsOrg2[geneyindex].split("\t")[1];
                                    
                                    if((index-indey)==1 || (index-indey)==-1){ //if the 2 genes are near each other in organism 2
                                        if(orix1.equalsIgnoreCase(oriy1)){
                                            if(!orix2.equalsIgnoreCase(oriy2)){
                                                noncontig++;
                                            }
                                            else if((index-indey)==-1 && orix1.equalsIgnoreCase(orix2)){
                                                contig++;
                                            }
                                            else if((index-indey)==1 && !orix1.equalsIgnoreCase(orix2)){
                                                contig++;
                                            }
                                            else{
                                                noncontig++;
                                            }
                                        }
                                        else{
                                            if(orix2.equalsIgnoreCase(oriy2)){
                                                noncontig++;
                                            }
                                            else if((index-indey)==-1 && orix1.equalsIgnoreCase(orix2)){
                                                contig++;
                                            }
                                            else if((index-indey)==1 && !orix1.equalsIgnoreCase(orix2)){
                                                contig++;
                                            }
                                            else{
                                                noncontig++;
                                            }
                                            
                                        }
                                    }
                                    else{
                                        noncontig++;
                                    }
                                }
                                int z = contig + noncontig;
                                double goc = (double) contig / reducedOrthologs.length;
                                DecimalFormat df = new DecimalFormat("#.###");
                                //System.out.println(contig + "\t" + noncontig + "\t" + reducedOrthologs.length + "\t" + df.format(goc));
                                collectGOC[u]=goc;

                            }
                            startnum+=increasenum;
                            //do statistics on goc values at a certain number of orthologs here (the values are in the array collectGOC)
                            mean(collectGOC);
                            stdev(collectGOC);
                        }
                    }
            }
            }
        }
    }
    public static void mean(double[] p) {
        double sum = 0;  // sum of all the elements
        for (int i = 0; i < p.length; i++) {
            sum += p[i];
        }
        DecimalFormat df = new DecimalFormat("#.###");
        double mean=sum/p.length;
        System.out.print(df.format(mean) + "\t");
    }
    //prints out the standard deviation of the sample
    public static void stdev(double[] p) {
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
        DecimalFormat df = new DecimalFormat("#.######");
        //double stdev=varsum/Math.sqrt(p.length-1); corrected stdev:
        double stdev=Math.sqrt(varsum/(p.length-1));
        System.out.println(df.format(stdev));
    }


}
