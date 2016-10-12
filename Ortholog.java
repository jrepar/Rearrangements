/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package reardist.GOCblastp.grimmgoctest;

/**
 *
 * @author Jelena
 */
public class Ortholog {
    private String pid1;
    private String pid2;
    private int pid1position;
    private int pid2position;
    private boolean plusStrand1;
    private boolean plusStrand2;
    private String ptt1;
    private String ptt2;
    private int name;
    
    private String gelement1;
    private String gelement2;
    private int gelementLength1;
    private int gelementLength2;
    private int geneNum1;
    private int geneNum2;

    public int getName() {
        return name;
    }

    public void setName(int name) {
        this.name = name;
    }
    
    public boolean isPlusStrand1() {
        return plusStrand1;
    }

    public void setPlusStrand1(boolean plusStrand) {
        this.plusStrand1 = plusStrand;
    }
    
    public boolean isPlusStrand2() {
        return plusStrand2;
    }

    public void setPlusStrand2(boolean plusStrand) {
        this.plusStrand2 = plusStrand;
    }
    
    public String getPid1() {
        return pid1;
    }

    public void setPid1(String pid1) {
        this.pid1 = pid1;
    }

    public String getPid2() {
        return pid2;
    }

    public void setPid2(String pid2) {
        this.pid2 = pid2;
    }

    public int getPid1position() {
        return pid1position;
    }

    public void setPid1position(int pid1position) {
        this.pid1position = pid1position;
    }

    public int getPid2position() {
        return pid2position;
    }

    public void setPid2position(int pid2position) {
        this.pid2position = pid2position;
    }

    public String getPtt1() {
        return ptt1;
    }

    public void setPtt1(String ptt1) {
        this.ptt1 = ptt1;
    }

    public String getPtt2() {
        return ptt2;
    }

    public void setPtt2(String ptt2) {
        this.ptt2 = ptt2;
    }

    public String getGelement1() {
        return gelement1;
    }

    public void setGelement1(String gelement1) {
        this.gelement1 = gelement1;
    }

    public String getGelement2() {
        return gelement2;
    }

    public void setGelement2(String gelement2) {
        this.gelement2 = gelement2;
    }

    public int getGelementLength1() {
        return gelementLength1;
    }

    public void setGelementLength1(int gelementLength1) {
        this.gelementLength1 = gelementLength1;
    }

    public int getGelementLength2() {
        return gelementLength2;
    }

    public void setGelementLength2(int gelementLength2) {
        this.gelementLength2 = gelementLength2;
    }

    public int getGeneNum1() {
        return geneNum1;
    }

    public void setGeneNum1(int geneNum1) {
        this.geneNum1 = geneNum1;
    }

    public int getGeneNum2() {
        return geneNum2;
    }

    public void setGeneNum2(int geneNum2) {
        this.geneNum2 = geneNum2;
    }
    
    
    
    
}
