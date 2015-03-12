package com.cndwzone.cluster;

import java.awt.Font;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.Scanner;
import java.util.Vector;

import org.jfree.chart.ChartColor;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.StandardChartTheme;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * @Description:Clustering by fast search and find of density peaks
 * @author drj
 * @date 2014-08-09
 * @version V1.0
 */

public class ClusterDemo {

	private static final double CLOCKS_PER_SEC = 1000;
	private static String datasetFilePath = "dataset1.txt";
	private String dirName; 
	private String outputHomePath="E:\\workspace\\Cluster\\"; 
	
	public static void main(String[] args) {

		long startTime, endTime;
		double executeTime;

		Vector<Point3D> pointsVector = new Vector<Point3D>();
		Vector<Point3D> descDeltaPointsVector = new Vector<Point3D>();
		Vector<Point3D> descRhoPointsVector = new Vector<Point3D>();
		Vector<Point3D> descRhoMultDeltaPointsVector = new Vector<Point3D>();
		Vector<Point3D> clusterCenterPointsVector = new Vector<Point3D>();
		
		Vector<Integer> rhoMultDelta = new Vector<Integer>();
		Vector<Integer> descRhoMultDelta = new Vector<Integer>();
		
		ArrayList<Double> allNumber = new ArrayList<Double>();
		
		ClusterDemo clusterDemo = new ClusterDemo();	
		clusterDemo.readFileAndChangeToNumber(datasetFilePath,allNumber);

		System.out.println("******start******");
		startTime=System.currentTimeMillis();   //获取开始时间
		
		clusterDemo.numberChangeToPoints(allNumber, pointsVector);
		double dc = clusterDemo.getDc(pointsVector, 0.016, 0.020);
		
		Vector<Integer> rho = clusterDemo.getLocalDensity(pointsVector, dc);
		Vector<Double> delta = clusterDemo.getNearestDistanceToHigherDensity(pointsVector, rho);
		Vector<Integer> tempRho = clusterDemo.getLocalDensity(pointsVector, dc);
		Vector<Double> tempDelta = clusterDemo.getNearestDistanceToHigherDensity(pointsVector, rho);
		
		clusterDemo.sortDescRhoAndSaveDescDeltaPointsVector(tempRho,descRhoPointsVector,pointsVector);
		clusterDemo.sortDescDeltaAndSaveDescDeltaPointsVector(tempDelta,descDeltaPointsVector,pointsVector);

		clusterDemo.paintShowDecisionTreeAndIdealFigure(delta,rho,pointsVector,rhoMultDelta,descRhoMultDelta, descRhoMultDeltaPointsVector);
		
		System.out.println("最佳dc值为：" + dc);
		System.out.println("nSamples:"+pointsVector.size());// 输出points的元素数
		
		System.out.println("看图分析，输入数据集可能的簇中心个数：");
		Scanner mayBeCenterNumber = new Scanner(System.in);
		int mayBeClusterCenterNum = mayBeCenterNumber.nextInt();
				
		Vector<Integer> descDeltaPointsVectorRho = clusterDemo.getLocalDensity(descDeltaPointsVector, dc);
		Vector<Double> descDeltaPointsVectorDelta = clusterDemo.getNearestDistanceToHigherDensity(descDeltaPointsVector, descDeltaPointsVectorRho);

		System.out.println("---通过delta确认的可能的簇中心---：");
		for (int i = 0; i < mayBeClusterCenterNum; i++) {
			String descDeltaPoints = clusterDemo.toString(descDeltaPointsVector.get(i));
			System.out.println("第"+(i+1)+"个可能的簇中心："+descDeltaPoints+"---descRho:"+descDeltaPointsVectorRho.get(i)+"---descDelta:"+descDeltaPointsVectorDelta.get(i));
		}
		System.out.println("---通过r确认的cluster point---：");
		for (int i = 0; i < mayBeClusterCenterNum; i++) {
			String descRhoAndDeltaPoints = clusterDemo.toString(descRhoMultDeltaPointsVector.get(i));
			System.out.println("第"+(i+1)+"个可能的簇中心："+descRhoAndDeltaPoints+" r="+descRhoMultDelta.get(i));
		}
		
		System.out.println("经过看图和输出的可能簇中心分析，输入数据集可能的簇中心个数：");
		Scanner centerNumber = new Scanner(System.in);
		int clusterCenterNum = centerNumber.nextInt();
		System.out.println("接下来，输入数据集可能的簇中心序列号：");
		ArrayList<Integer> clusterIndexArray = new ArrayList<Integer>();
		for (int i = 0; i < clusterCenterNum; i++) {
			Scanner centerIndex = new Scanner(System.in);
			int clusterIndex = centerIndex.nextInt();
			while (clusterIndex<1||clusterIndex>clusterCenterNum) {
				System.out.println("刚刚输入的序号不再其范围内，请重新输入：");
				clusterIndex = centerIndex.nextInt();
			}
			while (clusterIndexArray.contains(clusterIndex)) {				
				System.out.println("刚刚输入的序号已经存在，请重新输入：");
				clusterIndex = centerIndex.nextInt();
				while (clusterIndex<1||clusterIndex>clusterCenterNum) {
					System.out.println("刚刚输入的序号不再其范围内，请重新输入：");
					clusterIndex = centerIndex.nextInt();
				}
			}
			clusterIndexArray.add(clusterIndex);

		}
		
		clusterDemo.saveClusterCenterPointsVector(clusterCenterNum,clusterIndexArray,descRhoMultDeltaPointsVector,clusterCenterPointsVector);
		System.out.println("---通过看图和输出正式确认的簇中心---：");
		for (int i = 0; i < clusterCenterNum; i++) {
			String clusterCenterPoints = clusterDemo.toString(clusterCenterPointsVector.get(i));
			System.out.println("第"+(i+1)+"个簇中心："+clusterCenterPoints+" r="+descRhoMultDelta.get(i));
		}
		Vector<Point3D> rhoBigAndDistNearPointsVector = new Vector<Point3D>();
		Vector<Integer> descRhoPointsVectorRho = clusterDemo.getLocalDensity(descRhoPointsVector, dc);
		ArrayList<Integer> descRhoIndex = new ArrayList<Integer>();
		Vector<Double> descRhoPointsVectorDelta = clusterDemo.getNearestDistanceToHigherDensityPointsVector(descRhoPointsVector, descRhoPointsVectorRho, rhoBigAndDistNearPointsVector,descRhoIndex);

		@SuppressWarnings("unchecked")
		Vector<Point3D> clusterArray[] = new Vector[clusterCenterNum];
		Vector<Point3D> haloCluster = new Vector<Point3D>();
		
		clusterDemo.accordingToTheCenterForClusterCategories(clusterArray,clusterCenterNum,dc,clusterCenterPointsVector,descRhoPointsVector,descRhoPointsVectorDelta,descRhoPointsVectorRho,rhoBigAndDistNearPointsVector,descRhoIndex);
		
		clusterDemo.getHaloAndNewestCluster(clusterArray,haloCluster,dc);
		
		endTime=System.currentTimeMillis(); //获取结束时间
		executeTime =(endTime - startTime)/CLOCKS_PER_SEC;
		System.out.println("used time:"+executeTime);
		
	}
	
	/*
	 * point类型 x y z 转化成字符串类型(x,y,z)
	 */
	public String toString(Point3D point) {
		return "("+point.x+","+point.y+","+point.z+")";
	}

	/*
	 * 将txt文件中的数全部存放到ArrayList中
	 */
	private void readFileAndChangeToNumber(String datasetFilePath, ArrayList<Double> allNum) {
		try {
			String encoding = "GBK";
			File file = new File(datasetFilePath);
			if (file.isFile() && file.exists()) { // 判断文件是否存在
				InputStreamReader read = new InputStreamReader(new FileInputStream(file), encoding); // 考虑到编码格式
				BufferedReader bufferedReader = new BufferedReader(read);
				String lineTxt = null;
				while ((lineTxt = bufferedReader.readLine()) != null) {
					String lineTxtArray[] = lineTxt.split("\n");
					for (String oneLineTxt : lineTxtArray) {
						String oneLineTxtArray[] = oneLineTxt.split("	");
						for (String oneStringNum : oneLineTxtArray) {
							double oneNum = Double.parseDouble(oneStringNum);
							allNum.add(oneNum);
						}
					}
				}
				read.close();
			} else {
				System.out.println("Can not find the specified file!");
			}
		} catch (Exception e) {
			System.out.println("Error reading the file contents!");
			e.printStackTrace();
		}
		
	}

	/*
	 * 将ArrayList类型的allNumber中的数转成Point3D类型存入pointVector
	 */
	public void numberChangeToPoints(ArrayList<Double> allNumber ,Vector<Point3D> pointsVector ){
		
		for(int j=0;j<allNumber.size()-2;j++){	
			Point3D point = new Point3D(allNumber.get(j), allNumber.get(j+1), allNumber.get(j+2));
			j = j + 2;
			pointsVector.add(point);
	
		}
	}
	
	/*
	 * 求两点间距离
	 */
	private double getDistance(Point3D elementAt, Point3D elementAt2) {
		double tmp = Math.pow(elementAt.x - elementAt2.x, 2) + Math.pow(elementAt.y - elementAt2.y, 2) + Math.pow(elementAt.z - elementAt2.z, 2);
	    return Math.pow(tmp, 0.5);
	}
	
	/*
	 * 获取合适dc值
	 */
	private double getDc(Vector<Point3D> pointsVector, double neighborRateLow, double neighborRateHigh) {
		int nSamples = pointsVector.size();

		int nLow = (int) (neighborRateLow * nSamples * nSamples / 2);//这里没有理解，为什么是nSamples * nSamples-->现在已经明白
	    int nHigh = (int) (neighborRateHigh * nSamples * nSamples / 2);
	    double dc = 0.0;
	    int neighbors = 0;
	    System.out.println("nLow = "+nLow+", nHigh = "+nHigh);
	    while(neighbors < nLow || neighbors > nHigh){
	    	neighbors = 0;
	    	for(int i = 0; i < nSamples - 1; i++){// 两两计算距离
	            for(int j = i + 1; j < nSamples; j++){
	                if(getDistance(pointsVector.elementAt(i), pointsVector.elementAt(j)) <= dc)
	                	++neighbors;      // 如果距离小于等于dc   neighbors加1
	            }

	        }
	    	dc += 0.03;
	    	System.out.printf("dc = %.2f, neighbors = %d\n",dc, neighbors);
	    }
		return dc;
	}
	
	/*
	 * 局部密度
	 */
	private Vector<Integer> getLocalDensity(Vector<Point3D> pointsVector,double dc) {
		int nSamples = pointsVector.size();// 取data的元素数
	    Vector<Integer> rho = new Vector<Integer>();
	    for(int i = 0;i < nSamples;i++)//一个 int 类型n个元素，且值均为0的vecotr容器rho
	    	rho.add(0);
	    for(int i = 0; i < nSamples - 1; i++){
	        for(int j = i + 1; j < nSamples; j++){// 两两计算距离
	            if(getDistance(pointsVector.elementAt(i), pointsVector.elementAt(j)) < dc){
	            	int iNum = rho.get(i)+1;// 下标为i的元素值加1
	                rho.remove(i);
	                rho.add(i, iNum);

	                int jNum = rho.get(j)+1;// 下标为i的元素值加1
	                rho.remove(j);
	                rho.add(j, jNum);
	            }
	        }
	    }
	    return rho; 
	}
	
	/*
	 * 点到高局部密度点的最近距离
	 */
	private Vector<Double> getNearestDistanceToHigherDensity(Vector<Point3D> pointsVector, Vector<Integer> rho) {
		int nSamples = pointsVector.size();// 取points的元素数
		
		Vector<Double> delta = new Vector<Double>();
	    for(int i = 0; i < nSamples; i++){
	        double dist = 0.0;
	        boolean flag = false;
	        for(int j = 0; j < nSamples; j++){
	            if(i == j) continue;//如果是同一个点跳过
	            //如果j的局部密度比i大 然后求与i距离最小的距离
	            if(rho.get(j) > rho.get(i)){
	                double tmp = getDistance(pointsVector.get(i), pointsVector.get(j));
	                if(!flag){
	                    dist = tmp;
	                    flag = true;
	                }else {
	                	dist = tmp < dist ? tmp : dist; //dist记录两个点之间最小值

	                }
	            }
	        }	        
	        //如果j都比i的局部密度小，说明i的局部密度最大，则dist为i的是与j最大的距离
	        if(!flag){
	            for(int j = 0; j < nSamples; j++){
	                double tmp = getDistance(pointsVector.get(i), pointsVector.get(j));
	                dist = tmp > dist ? tmp : dist;
	            }
	        }
	        delta.add(i, dist);
	    }
	    return delta;
	}
	
	/*
	 * 点到高局部密度点的最近距离，并保存相应的坐标点和对应点的标号
	 */
	private Vector<Double> getNearestDistanceToHigherDensityPointsVector(Vector<Point3D> pointsVector, 
			Vector<Integer> rho,Vector<Point3D> rhoBigAndDistNearPointsVector,ArrayList<Integer> descRhoIndex) {
		int nSamples = pointsVector.size();// 取points的元素数
		
		Vector<Double> delta = new Vector<Double>();
	    for(int i = 0; i < nSamples; i++){
	        double dist = 0.0;
	        int tmpIndex = nSamples-1;
	        boolean flag = false;
	        for(int j = 0; j < nSamples; j++){
	            if(i == j) continue;//如果是同一个点跳过
	            //如果j的局部密度比i大 然后求与i距离最小的距离
	            if(rho.get(j) > rho.get(i)){
	                double tmp = getDistance(pointsVector.get(i), pointsVector.get(j));
	                if(!flag){
	                    dist = tmp;
	                    tmpIndex = j;
	                    flag = true;
	                }else {
	                	if(tmp < dist){
	                		dist = tmp;
	                		tmpIndex = j;
	                	}
	                }
	            }
	        }
	        rhoBigAndDistNearPointsVector.add(pointsVector.get(tmpIndex));
	        descRhoIndex.add(tmpIndex);
	        
	        //如果j都比i的局部密度小，说明i的局部密度最大，则dist为i的是与j最大的距离
	        if(!flag){
	            for(int j = 0; j < nSamples; j++){
	                double tmp = getDistance(pointsVector.get(i), pointsVector.get(j));
	                dist = tmp > dist ? tmp : dist;
	            }
	        }
	        delta.add(i, dist);
	    }
	    return delta;
	}
	
	/*
	 * 对rho按照降序排序，并且生成新的对应坐标点
	 */
	private void sortDescRhoAndSaveDescDeltaPointsVector(Vector<Integer> tempRho, Vector<Point3D> descDeltaPointsVector, Vector<Point3D> pointsVector) {
		//对delta的值sort
		Vector<Integer> descRho = new Vector<Integer>();
		for(int i = 0;i < tempRho.size();i++){
			double temp = tempRho.get(i);
			int tempIndex = i;
			for(int j = 0;j < tempRho.size();j++){
				if(temp<tempRho.get(j)){
					temp = tempRho.get(j);
					tempIndex = j;
				}
			}
			descRho.add(tempRho.get(tempIndex));
			tempRho.remove(tempIndex);
			tempRho.add(tempIndex, 0);
			descDeltaPointsVector.add(pointsVector.get(tempIndex));
		}
		
	}

	/*
	 * 对delta按照降序排序，并且生成新的对应坐标点
	 */
	private void sortDescDeltaAndSaveDescDeltaPointsVector(Vector<Double> tempDelta, Vector<Point3D> descDeltaPointsVector, Vector<Point3D> pointsVector) {
		//对delta的值sort
		Vector<Double> descDelta = new Vector<Double>();
		for(int i = 0;i < tempDelta.size();i++){
			double temp = tempDelta.get(i);
			int tempIndex = i;
			for(int j = 0;j < tempDelta.size();j++){
				if(temp<tempDelta.get(j)){
					temp = tempDelta.get(j);
					tempIndex = j;
				}
			}
			descDelta.add(tempDelta.get(tempIndex));
			tempDelta.remove(tempIndex);
			tempDelta.add(tempIndex, 0.0);
			descDeltaPointsVector.add(pointsVector.get(tempIndex));
		}
		
	}
	
	/*
	 * 对rho和delta的积降序排序
	 */
	private void descRhoMultDelta(Vector<Integer> rhoMultDelta, Vector<Integer> descRhoMultDelta, Vector<Point3D> descRhoMultDeltaPointsVector, Vector<Point3D> pointsVector) {
		for(int i = 0;i < rhoMultDelta.size();i++){
			double temp = rhoMultDelta.get(i);
			int tempIndex = i;
			for(int j = 0;j < rhoMultDelta.size();j++){
				if(temp<rhoMultDelta.get(j)){
					temp = rhoMultDelta.get(j);
					tempIndex = j;
				}
			}
			descRhoMultDelta.add(rhoMultDelta.get(tempIndex));
			rhoMultDelta.remove(tempIndex);
			rhoMultDelta.add(tempIndex, 0);
			descRhoMultDeltaPointsVector.add(pointsVector.get(tempIndex));
		}
		
	}
	
	/*
	 * 根据rho*delta值的大小人为判断簇中心个数，保存对应簇中心坐标
	 */
//	private void mayBeClusterCenterPointsVector(int mayBeClusterCenterNum,Vector<Point3D> descRhoMultDeltaPointsVector, Vector<Point3D> mayBeClusterCenterPointsVector) {
//		for (int i = 0; i < mayBeClusterCenterNum; i++) {
//			mayBeClusterCenterPointsVector.add(descRhoMultDeltaPointsVector.get(i));
//		}
//		
//	}
	
	/*
	 * 保存选择的簇中心
	 */
	private void saveClusterCenterPointsVector(int clusterCenterNum,
			ArrayList<Integer> clusterIndexArray,
			Vector<Point3D> descRhoMultDeltaPointsVector,
			Vector<Point3D> clusterCenterPointsVector) {
		
		for (int i = 0; i < clusterCenterNum; i++) {
			clusterCenterPointsVector.add(descRhoMultDeltaPointsVector.get(clusterIndexArray.get(i)-1));
		}
	}

	/*
	 * 判断i不是簇中心点
	 */
	private boolean contains(ArrayList<Point3D> arrayListCenter, Point3D point3d) {
		
		for (Point3D vector : arrayListCenter) {
			 if ( vector.equals(point3d) )
				 return true;	
		}
	    return false;
	}
	
	/*
	 * 物以类聚
	 */
	private void accordingToTheCenterForClusterCategories(Vector<Point3D>[] clusterArray,
			int clusterCenterNum, 
			double dc, Vector<Point3D> clusterCenterPointsVector,
			Vector<Point3D> descRhoPointsVector, 
			Vector<Double> descRhoPointsVectorDelta, Vector<Integer> descRhoPointsVectorRho,
			Vector<Point3D> rhoBigAndDistNearPointsVector, ArrayList<Integer> descRhoIndex) {
		
		int nSamples = descRhoPointsVector.size();// 取points的元素数
		for(int i=0;i<clusterCenterNum;i++){
			clusterArray[i]=new Vector<Point3D>();
			clusterArray[i].add(clusterCenterPointsVector.get(i));
		}
		
		ArrayList<Point3D> arrayListCenter = new ArrayList<Point3D>();
		for(int i=0;i<clusterCenterNum;i++){
			arrayListCenter.add(clusterCenterPointsVector.get(i));
		}
		
		System.out.println("nSamples:"+nSamples);
	    for(int i = nSamples-1; i >= 0; i--){
	    	ArrayList<Integer> roadIndexArray = new ArrayList<Integer>();
	    	int temp = i;
	    	while(!contains(arrayListCenter,descRhoPointsVector.get(temp))){
	    		for (int j = 0; j < clusterCenterNum; j++) {
	    			if(rhoBigAndDistNearPointsVector.get(temp) == clusterArray[j].get(0)){
	    				String descRhoPoint = toString(descRhoPointsVector.get(i));
	    				String rhoBigAndDistNearPoint = toString(rhoBigAndDistNearPointsVector.get(i));
	    				System.out.println("簇"+(j+1)+"添加的是第:"+(i+1)+"点对应的坐标是："+descRhoPoint+"对应密度大距离近点的坐标是："+rhoBigAndDistNearPoint);
	    				
		    			clusterArray[j].add(descRhoPointsVector.get(i));
		    			
	    				System.out.print("寻找点对应簇中心路线显示:");
	    				System.out.print(descRhoPoint+"-->");
	    				for (Integer roadIndex : roadIndexArray) {
	    					String findClusterCenterRoadPoint = toString(descRhoPointsVector.get(roadIndex));
							System.out.print(findClusterCenterRoadPoint+"-->");
						}
	    				
	    				String clusterCenterPoint = toString(rhoBigAndDistNearPointsVector.get(temp));
	    				System.out.println(clusterCenterPoint);

	    			}
				}
	    		temp = descRhoIndex.get(temp);
	    		roadIndexArray.add(temp);
		    }

	    }
//	    paintShowCluster(clusterArray);
	}
	
	/*
	 * 将halo光晕分离出来
	 */
	private void getHaloAndNewestCluster(Vector<Point3D>[] clusterArray,
			Vector<Point3D> haloCluster, double dc) {

		@SuppressWarnings("unchecked")
		Vector<Double> borderRho[]=new Vector[clusterArray.length];
		for(int i=0;i<clusterArray.length;i++){
			borderRho[i]=new Vector<Double>();
			borderRho[i].add(0.0);
		}
		
		if (clusterArray.length > 1) {
			for (int i = 0; i < clusterArray.length-1; i++) {
				Vector<Integer> rhoFirst = getLocalDensity(clusterArray[i],dc);
				for (int i1 = 0; i1 < clusterArray[i].size(); i1++) {
					for (int j = i+1; j < clusterArray.length; j++) {
						Vector<Integer> rhoSecond = getLocalDensity(clusterArray[j],dc);
						for (int j2 = 0; j2 < clusterArray[j].size(); j2++) {
							double dist = getDistance(clusterArray[i].get(i1),clusterArray[j].get(j2) );
							if(dist <= dc){
								double rhoAver=(rhoFirst.get(i1)+rhoSecond.get(j2))/2;
								if (rhoAver > borderRho[i].get(0)) {
									borderRho[i].remove(0);
									borderRho[i].add(0, rhoAver);
								}
								if (rhoAver > borderRho[j].get(0)) {
									borderRho[j].remove(0);
									borderRho[j].add(0, rhoAver);
								}
								String clusterArrayIpoint = toString(clusterArray[i].get(i1));
								String clusterArrayJpoint = toString(clusterArray[j].get(j2));
								System.out.println("不同簇但两点距离小于dc的点的坐标："+clusterArrayIpoint+"和"+clusterArrayJpoint);
							}
						}
					}
				}
				
			}
			int i1 = 0;
			int j1 = 0;
			for (int i = 0; i < clusterArray.length; i++) {
				Vector<Integer> rho = getLocalDensity(clusterArray[i],dc);
				
				for (int j = 0; j < clusterArray[i].size(); j++) {
					if (rho.get(j) < borderRho[i].get(0)) {
						haloCluster.add(clusterArray[i].get(j));
						clusterArray[i].remove(j);
						if(i == 0){
							i1++;
						}
						if(i == 1){
							j1++;
						}
					}

				}
				
			}
			ArrayList<Integer> haloNum = new ArrayList<Integer>();
			haloNum.add(i1);
			haloNum.add(j1);
			
			for (int i = 0; i < borderRho.length; i++) {
				System.out.println("簇"+(i+1)+"的halo局部密度:"+borderRho[i].get(0)+" 光晕数："+haloNum.get(i));
			}
		}
		
		paintShowCluster(clusterArray,haloCluster);
	}	
	
	/*
	 * 物以类聚结果图
	 */
	private void paintShowCluster(Vector<Point3D>[] nonHaloClusterArray,Vector<Point3D> haloCluster) {

		XYSeries xyseries1 = new XYSeries("簇1坐标点");
		XYSeries xyseries2 = new XYSeries("簇2坐标点");
		XYSeries xyseries3 = new XYSeries("光晕坐标点");
	
			for (int i = 0; i < nonHaloClusterArray[0].size(); i++) {
				xyseries1.add(nonHaloClusterArray[0].get(i).x, nonHaloClusterArray[0].get(i).y);
			}
			
			for (int i = 0; i < nonHaloClusterArray[1].size(); i++) {
				xyseries2.add(nonHaloClusterArray[1].get(i).x, nonHaloClusterArray[1].get(i).y);
			}
			
			for (int i = 0; i < haloCluster.size(); i++) {
				xyseries3.add(haloCluster.get(i).x, haloCluster.get(i).y);
			}	
				
		XYSeriesCollection xyseriescollection1 = new XYSeriesCollection(); 
		xyseriescollection1.addSeries(xyseries1);
		xyseriescollection1.addSeries(xyseries2);
		xyseriescollection1.addSeries(xyseries3);

		StandardChartTheme standardChartTheme = new StandardChartTheme("CN");
		standardChartTheme.setExtraLargeFont(new Font("隶书", Font.BOLD, 20));
		standardChartTheme.setRegularFont(new Font("宋书", Font.PLAIN, 15));
		standardChartTheme.setLargeFont(new Font("宋书", Font.PLAIN, 15));
		ChartFactory.setChartTheme(standardChartTheme);

		JFreeChart chart1 = ChartFactory.createScatterPlot("点分布图", "X轴","Y轴", xyseriescollection1, PlotOrientation.VERTICAL, true,false, false);

		XYPlot p = chart1.getXYPlot();
		// 设置图的背景颜色
		p.setBackgroundPaint(ChartColor.WHITE);			
		
		try {
			ChartUtilities.saveChartAsPNG(new File(outputHomePath + dirName+"\\real分布图.png"), chart1, 1000, 1000);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	

	/*
	 * 将局部密度和距离显示在决策图上，并且生成对应的r值图
	 */
	private void paintShowDecisionTreeAndIdealFigure(Vector<Double> delta, Vector<Integer> rho,Vector<Point3D> pointsVector, Vector<Integer> rhoMultDelta ,Vector<Integer> descRhoMultDelta, Vector<Point3D> descRhoMultDeltaPointsVector) {

		XYSeries xyseries = new XYSeries("局部密度和到高局部密度点的距离");
		for (int i = 0; i < delta.size() ; i++) {
			xyseries.add(rho.get(i), delta.get(i));
			rhoMultDelta.add((int) (rho.get(i)*delta.get(i)));
		}
		
		descRhoMultDelta(rhoMultDelta,descRhoMultDelta,descRhoMultDeltaPointsVector,pointsVector);
		
		XYSeries xyseries1 = new XYSeries("簇1坐标点");
		XYSeries xyseries2 = new XYSeries("簇2坐标点");
		XYSeries xyseries3 = new XYSeries("簇3坐标点");
		XYSeries xyseries4 = new XYSeries("簇4坐标点");
		XYSeries xyseries5 = new XYSeries("簇5坐标点");
		XYSeries xyseries6 = new XYSeries("簇6坐标点");
		XYSeries xyseries7 = new XYSeries("簇7坐标点");
		
		for (int j = 0; j < pointsVector.size() ; j++) {
			if(pointsVector.get(j).z == 1){
				xyseries1.add(pointsVector.get(j).x, pointsVector.get(j).y);
			}else if(pointsVector.get(j).z == 2){
				xyseries2.add(pointsVector.get(j).x, pointsVector.get(j).y);
			}else if(pointsVector.get(j).z == 3){
				xyseries3.add(pointsVector.get(j).x, pointsVector.get(j).y);
			}else if(pointsVector.get(j).z == 4){
				xyseries4.add(pointsVector.get(j).x, pointsVector.get(j).y);
			}else if(pointsVector.get(j).z == 5){
				xyseries5.add(pointsVector.get(j).x, pointsVector.get(j).y);
			}else if(pointsVector.get(j).z == 6){
				xyseries6.add(pointsVector.get(j).x, pointsVector.get(j).y);
			}else{
				xyseries7.add(pointsVector.get(j).x, pointsVector.get(j).y);
			}
		}
		
		XYSeries xyseries8 = new XYSeries("r值");
		for (int i = 0; i < descRhoMultDelta.size() ; i++) {
			xyseries8.add(i, descRhoMultDelta.get(i));
		}
			
		XYSeriesCollection xyseriescollection = new XYSeriesCollection(); // 再用XYSeriesCollection添加入XYSeries对象
		xyseriescollection.addSeries(xyseries);
		
		XYSeriesCollection xyseriescollection1 = new XYSeriesCollection(); 
		xyseriescollection1.addSeries(xyseries1);
		xyseriescollection1.addSeries(xyseries2);
		xyseriescollection1.addSeries(xyseries3);
		xyseriescollection1.addSeries(xyseries4);
		xyseriescollection1.addSeries(xyseries5);
		xyseriescollection1.addSeries(xyseries6);
		xyseriescollection1.addSeries(xyseries7);

		XYSeriesCollection xyseriescollection2 = new XYSeriesCollection(); 
		xyseriescollection2.addSeries(xyseries8);
		
		StandardChartTheme standardChartTheme = new StandardChartTheme("CN");
		standardChartTheme.setExtraLargeFont(new Font("隶书", Font.BOLD, 20));
		standardChartTheme.setRegularFont(new Font("宋书", Font.PLAIN, 15));
		standardChartTheme.setLargeFont(new Font("宋书", Font.PLAIN, 15));
		ChartFactory.setChartTheme(standardChartTheme);

		JFreeChart chart = ChartFactory.createScatterPlot("决策图(decision tree)", "ρ轴","δ轴", xyseriescollection, PlotOrientation.VERTICAL, true,false, false);
		JFreeChart chart1 = ChartFactory.createScatterPlot("点分布图", "X轴","Y轴", xyseriescollection1, PlotOrientation.VERTICAL, true,false, false);
		JFreeChart chart2 = ChartFactory.createScatterPlot("r分布图", "n","r", xyseriescollection2, PlotOrientation.VERTICAL, true,false, false);

		XYPlot p = chart1.getXYPlot();
		
		p.setBackgroundPaint(ChartColor.WHITE); // 设置图的背景颜色
			
		SimpleDateFormat df = new SimpleDateFormat("yyyyMMddHHmmss");
		
		dirName = "img\\"+df.format(new Date());  

		createDir(dirName);
		
		try {
			ChartUtilities.saveChartAsPNG(new File(outputHomePath+dirName+"\\决策图.png"), chart, 1000, 1000);
			ChartUtilities.saveChartAsPNG(new File(outputHomePath+dirName+"\\理想分布图.png"), chart1, 1000, 1000);
			ChartUtilities.saveChartAsPNG(new File(outputHomePath+dirName+"\\r值图.png"), chart2, 10000, 1000);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/*
	 * 创建图片存放文件夹
	 */
	public static boolean createDir(String destDirName) {  
	    File dir = new File(destDirName);  
	    if(dir.exists()) {  
		    System.out.println("Create a directory: " + destDirName + " Failed, the target directory already exists!");  
		    return false;  
	    }  
	    if(dir.mkdirs()) {  
	    	System.out.println("Create a directory: " + destDirName + " successfully!");  
	    	return true;  
	    } else {  
	    	System.out.println("Create a directory: " + destDirName + " Failed!");  
	    	return false;  
	    }  
	} 
		
}

class Point3D {
    double x;
    double y;
    double z;
    
	public Point3D(double x, double y, double z) {
		super();
		this.x = x;
		this.y = y;
		this.z = z;
	}

}
