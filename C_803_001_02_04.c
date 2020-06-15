#include <stdio.h>
#include <math.h>
#include <float.h>


#define pi 3.141592653589
#define length 840	//总长度840mm，绝热段100mm，蒸发段冷凝段每段110mm

/*函数原型*/
double mu_l_given_T(double T);
double cp_l_given_T(double T);
double cv_l_given_T(double T);
double lamt_l_given_T(double T);
double Pr_l_given_T(double T);
double ro_l_given_T(double T);
double p_sat_given_T(double T);
double T_sat_given_p(double p);
double ro_sat_given_T(double T);
double h_fg_given_T(double T);
double u_given_vt(double volume, double mass, double temp);
double T_given_uv(double u, double volume, double mass, double T1);


int main()
{
	
	/* 定义步长 */
	double dx=1e-3;      /*距离步长，m */
	double dt=2e-7;      /*时间步长,s*/

	/* 定义总时间-输出步长间隔 */
	double T_total=20;	//总时长，s
	long nn=T_total/dt;	//总时间节点数
	int print_fre=5000;	//输出的步长间隔，每5000个时间步长输出一次（每0.001s输出一次）
	long heat_print_fre=0.001/dt;     /*输出换热量的步长间隔（毎0.001s输出一次平均换热量）*/

	/* 几何参数 */
	const double d = 2.3e-3;	//PHP内径，m
	const double d_out = 3.3e-3;	//PHP外径，m
	const double film = 5.2e-5;	//液膜厚度，m
	const double A_v = 0.25*pi*pow((d - 2 * film), 2); //气泡截面积，m2
	const double A_l=4.15476e-6;	//液塞截面积，m2，内径2.3mm的圆
	const double A_w = 4.39823e-6;	//管壁截面积，m2，内径2.3mm外径3.3mm的圆环

	const double g=9.8;	//重力加速度，m/s2
	
	/* 物性参数 */
	const double ro_l = 72.26;	//液氢19K密度，kg/m3
	const double mm = 4124;    //氢气的气体常数EES，J/kg-K
	const double lamt_w = 2.169;                           /*20K,SS304热导率，W/m-K*/
	const double c_w = 13.45;                            /*20K,SS304热容,J/kg-K*/
	const double ro_w = 8072;                           /*20K,SS304管壁的密度，kg/m3*/

	const double T_min=18;	//温度下限，K
	const double T_max=31;	//温度上限，K
	/* 管外温度（预留量）*/
	const double T_out = 20;
	const double T_c = 19.065;	//固定冷凝端温度，K
	
	/* 状态值 */
	int codenumber = 0;	//0表示正常，1表示（n=1气体开始）液塞移动过快或时间步长过大，2表示（n=1液体开始）液塞移动过快或时间不长过大
	long bubble_gen_total = 0;	//新产生气泡的总数量


    /* ****************定义控制体的诸多变量*******************/
	int sort[length+1];	//控制体分类数组，数组下标为0-length，所以下标长度为length+1，=1表示液，=2表示气，=3表示液-气，=4表示气-液，=8表示液-气-液，=9表示液塞要合并，=6表示新气泡产生
	int c_gravity[length+1];	//重力方向系数，值为-1，0，1
	int plain[length+1];	//加热形式标记，1表示设定壁温，0表示壁温可变
	float heat_flux[length+1];	//热流密度，W/m
	int heat_sig[length+1];		//传热标记，1表示蒸发段，0表示绝热段，-1表示冷凝段
	float h_v[length+1];	//气泡中的蒸发冷凝换热系数，W/m2-K
	float h_c[length+1];	//管外对流换热系数？设为0
	double h_b[length+1];	//沸腾换热系数，W/m2-K
	int Tr_v_sig[length+1];	//液膜相变标记，1蒸发，0无相变，-1冷凝
	double h_lleft[length+1];	//液体单相对流换热系数，W/m2-K
	double h_lright[length+1];	//液体单相对流换热系数，W/m2-K

	double mu_left[length+1];	//每个控制体内左边液塞动力粘度，Pa-s
	double mu_right[length+1];	//每个控制体内右边液塞动力粘度，Pa-s
	double c_pl_left[length+1];	//每个控制体内左边液塞比热，J/kg-K
	double c_pl_right[length+1];	//每个控制体内右边液塞比热，J/kg-K
	double lamt_left[length+1];	//每个控制体内左边液塞热导率，W/m-K
	double lamt_right[length+1];	//每个控制体内右边液塞热导率，W/m-K
	double h_fg_lleft[length+1];	//每个控制体内左边液塞汽化潜热，J/kg
	double h_fg_lright[length+1];	//每个控制体内右边液塞汽化潜热，J/kg
	double h_fg_v[length+1];	//控制体内气体部分的汽化潜热，与液体的汽化潜热相等
	double C_l[length+1];	//控制体的液体阻力系数

	double x_lleft[length+1][3];	//每个控制体内左边液塞长度，m
	double x_lright[length+1][3];	//每个控制体内右边液塞长度，m
	double v_lleft[length+1][3];	//每个控制体内左边液塞速度，m/s
	double v_lright[length+1][3];	//每个控制体内右边液塞速度，m/s
	double T_lleft[length+1][3];	//每个控制体内左边液塞温度，m/s
	double T_lright[length+1][3];	//每个控制体内右边液塞温度，m/s
	double P_l[length+1][3];	//液塞压力，Pa
	double x_v[length+1][3];	//控制体中气塞长度，m
	double ro_v[length+1][3];	//控制体中气塞密度，kg/m3
	double T_v[length+1][3];	//控制体中气塞温度，K
	double T_sat[length + 1][3];	//控制体中液膜表面饱和温度，K
	double P_v[length+1][3];	//控制体中气塞压力，Pa

	double T_w[length+1][3];	//壁面温度，K

	/* 初始化气液分布，10个气塞10个液塞，充液率50%*/
	int n;	//n表示控制体序号（位置坐标），从1到840，每个控制体长1mm
	{
		sort[0] = 0;	//数组的第一个元素用不到，赋值0
		for (n = 1; n <= 20; n = n + 1)
		{sort[n] = 2;}	//气塞控制体
		sort[21] = 4;	//气-液控制体
		for (n = 22; n <= 62; n = n + 1)
		{sort[n] = 1;}	//液塞控制体
		sort[63] = 3;	//液-气控制体

		for (n = 64; n <= 104; n = n + 1)
		{sort[n] = 2;}
		sort[105] = 4;
		for (n = 106; n <= 146; n = n + 1)
		{sort[n] = 1;}
		sort[147] = 3;

		for (n = 148; n <= 188; n = n + 1)
		{sort[n] = 2;}
		sort[189] = 4;
		for (n = 190; n <= 230; n = n + 1)
		{sort[n] = 1;}
		sort[231] = 3;

		for (n = 232; n <= 272; n = n + 1)
		{sort[n] = 2;}
		sort[273] = 4;
		for (n = 274; n <= 314; n = n + 1)
		{sort[n] = 1;}
		sort[315] = 3;

		for (n = 316; n <= 356; n = n + 1)
		{sort[n] = 2;}
		sort[357] = 4;
		for (n = 358; n <= 398; n = n + 1)
		{sort[n] = 1;}
		sort[399] = 3;

		for (n = 400; n <= 440; n = n + 1)
		{sort[n] = 2;}
		sort[441] = 4;
		for (n = 442; n <= 482; n = n + 1)
		{sort[n] = 1;}
		sort[483] = 3;

		for (n = 484; n <= 524; n = n + 1)
		{sort[n] = 2;}
		sort[525] = 4;
		for (n = 526; n <= 566; n = n + 1)
		{sort[n] = 1;}
		sort[567] = 3;

		for (n = 568; n <= 608; n = n + 1)
		{sort[n] = 2;}
		sort[609] = 4;
		for (n = 610; n <= 650; n = n + 1)
		{sort[n] = 1;}
		sort[651] = 3;

		for (n = 652; n <= 692; n = n + 1)
		{sort[n] = 2;}
		sort[693] = 4;
		for (n = 694; n <= 734; n = n + 1)
		{sort[n] = 1;}
		sort[735] = 3;

		for (n = 736; n <= 776; n = n + 1)
		{sort[n] = 2;}
		sort[777] = 4;
		for (n = 778; n <= 818; n = n + 1)
		{sort[n] = 1;}
		sort[819] = 3;
		for (n = 820; n <= 840; n = n + 1)
		{sort[n] = 2;}
	}
	/* 设定重力方向系数*/
	{
		for (n = 1; n <= 155; n = n + 1)
		{c_gravity[n] = 1;}
		for (n = 156; n <= 365; n = n + 1)
		{c_gravity[n] = -1;}
		for (n = 366; n <= 575; n = n + 1)
		{c_gravity[n] = 1;}
		for (n = 576; n <= 785; n = n + 1)
		{c_gravity[n] = -1;}
		for (n = 786; n <= 840; n = n + 1)
		{c_gravity[n] = 1;}
	}
	/* 设定壁温固定系数*/
	{
		for (n = 1; n <= length; n = n + 1)
		{
			plain[n] = 0;	//壁温可变
		}
		for (n = 311; n <= 420; n = n + 1)
		{
			plain[n] = 1;	//壁温被设定，冷凝段1
		}
		for (n = 731; n <= 840; n = n + 1)
		{
			plain[n] = 1;	//壁温被设定，冷凝段2
		}
	}
	/* 设定热流量(W/m)和蒸发段1和2的位置*/
	float heat_flux_set = 0.4 / 0.11;	//蒸发段热流量，0.8W平均分配到两段蒸发段共0.11m的长度上
	int n_heat_start1 = 101;	//蒸发段起始控制体序号
	int n_heat_end1 = 210;	//蒸发段终止控制体序号
	int n_heat_start2 = 521;	//蒸发段起始控制体序号
	int n_heat_end2 = 630;	//蒸发段终止控制体序号
	/*设定各部分控制体的热流密度和传热标记*/
	{
		for (n = 1; n <= 100; n = n + 1)	//绝热段1
		{	heat_flux[n] = 0;
			heat_sig[n] = 0;		}
		for (n = 101; n <= 210; n = n + 1)	//蒸发段1
		{	heat_flux[n] = heat_flux_set;
			heat_sig[n] = 1;		}
		for (n = 211; n <= 310; n = n + 1)	//绝热段2
		{	heat_flux[n] = 0;
			heat_sig[n] = 0;		}
		for (n = 311; n <= 420; n = n + 1)	//冷凝段1，热流密度为0
		{	heat_flux[n] = 0;
			heat_sig[n] = -1;		}
		for (n = 421; n <= 520; n = n + 1)	//绝热段3
		{	heat_flux[n] = 0;
			heat_sig[n] = 0;		}
		for (n = 521; n <= 630; n = n + 1)	//蒸发段2
		{	heat_flux[n] = heat_flux_set;
			heat_sig[n] = 1;		}
		for (n = 631; n <= 730; n = n + 1)	//绝热段4
		{	heat_flux[n] = 0;
			heat_sig[n] = 0;		}
		for (n = 731; n <= 840; n = n + 1)	//冷凝段2
		{	heat_flux[n] = 0;
			heat_sig[n] = -1;		}
	}
	/* 设定气泡过热度，K*/
	float overheat = 0.25;
	/* 设定气泡产生时间间隔*/
	double generate_fre = 0.05;	//气泡产生的时间间隔，s
	double generate_inter1 = 0;	//第1个气泡产生点的计时器
	double generate_inter2 = 0;	//第2个气泡产生点的计时器
	double generate_inter3 = 0;	//第3个气泡产生点的计时器
	double generate_inter4 = 0;	//第4个气泡产生点的计时器
	double generate_inter5 = 0;	//第5个气泡产生点的计时器
	double generate_inter6 = 0;	//第6个气泡产生点的计时器
	double generate_inter7 = 0;	//第7个气泡产生点的计时器
	double generate_inter8 = 0;	//第8个气泡产生点的计时器
	/* 设定产生、消失最小长度*/
	double l_disappear = 0.2*dx;
	/* 设定气泡中蒸发与冷凝换热系数*/
	for (n = 1; n <= length; n = n + 1)
	{
		h_v[n] = 0;
		if (heat_sig[n] == 1)
		{h_v[n] = 2028;}	//蒸发换热系数1000W/m2-K
		if (heat_sig[n] == -1)
		{h_v[n] = 2028;}	//冷凝换热系数1000W/m2-K
	}
	/* 设定管外冷却换热参数（预留）*/
	for (n = 1; n <= length; n = n + 1)
	{h_c[n] = 0;}
	/* 设定初始相变标记*/
	for (n = 1; n <= length; n = n + 1)
	{Tr_v_sig[n] = 0;	//初始设定为无相变
	}
	/*计算交界面个数及位置*/
	int transfer_point_number;	//用来记录气液界面序号,累加变量
	int	n_trans_total;	//气液界面总数
	int n_bubble_total, n_liquid_total;	//用来记录气泡、液塞总数
	int n_trans[100];	//用来记录气液界面位置坐标
	{
		n = 1;
		if (sort[n] == 3 || sort[n] == 4)	//液气或气液控制体
		{
			transfer_point_number = 1;
			n_trans[1] = 1;	//第[1]个气液界面的位置为1
		}
		else
		{
			transfer_point_number = 0;
		}              /*表示气泡和液塞交界面的数量*/
		n = 2;
		while (n <= length)
		{
			if (sort[n] == 3 || sort[n] == 4)
			{
				transfer_point_number = transfer_point_number + 1;
				n_trans[transfer_point_number] = n;                  /*使用n_trans（）来记录每个交界面的位置*/
			}
			n = n + 1;
		}
		n_trans_total = transfer_point_number;                 /*记录交界面总数*/
		n_bubble_total = n_trans_total / 2;                      /*气泡的总数*/
		n_liquid_total = n_trans_total / 2;                      /*液塞的总数*/
	}
	
    /* 初始化控制体数据part01（气塞和液塞的长度）*/
	n=1;
	while (n<=length)
	{
		if (sort[n]==3)	//液-气控制体，左边一半是液，右边一半是气 
		{
			x_lleft[n][1]=dx/2;	//数组下标n表示控制体序号，1表示时间步数
			x_lright[n][1]=0;
			x_v[n][1]=dx/2;
		}
		if (sort[n]==4)	//气-液控制体，左边一半是气，右边一半是液
		{
			x_lleft[n][1]=0;
			x_lright[n][1]=dx/2;
			x_v[n][1]=dx/2;
		}
		if (sort[n]==1)	//液控制体
		{
			x_lleft[n][1]=dx;	//？？？为啥不是dx/2
			x_lright[n][1]=dx;	//???为啥不是dx/2
			x_v[n][1]=0;
		}
		if (sort[n]==2)	//气控制体
		{
			x_lleft[n][1]=0;
			x_lright[n][1]=0;
			x_v[n][1]=dx;
		}
		n=n+1;
	}
	/* 初始化控制体数据part02（液体速度，温度，粘度，比热，热导率，汽化潜热，
	气体的密度，温度，压力，汽化潜热,壁面温度，液体压力）*/
	for (n=1;n<=length;n=n+1)
	{
		if (sort[n]==1)	//液控制体
		{
			v_lleft[n][1] = 0;	//初始速度为0
			v_lright[n][1] = 0;
			T_lleft[n][1] = T_c;	//初始温度与冷凝段温度相同
			T_lright[n][1] = T_c;

			mu_left[n]=mu_l_given_T(T_lleft[n][1]);
			mu_right[n]=mu_l_given_T(T_lright[n][1]);
			c_pl_left[n]=cp_l_given_T(T_lleft[n][1]);
			c_pl_right[n]=cp_l_given_T(T_lright[n][1]);
			lamt_left[n]=lamt_l_given_T(T_lleft[n][1]);
			lamt_right[n]=lamt_l_given_T(T_lright[n][1]);
			h_fg_lleft[n]=h_fg_given_T(T_lleft[n][1]);
			h_fg_lright[n]= h_fg_given_T(T_lright[n][1]);
			
			ro_v[n][1]=0;
			T_v[n][1]=0;
			T_sat[n][1] = 0;
			P_v[n][1]=0;
			h_fg_v[n]=0;
		}

		if (sort[n]==2)	//气控制体
		{
			v_lleft[n][1] = 0;
			v_lright[n][1] = 0;
			T_lleft[n][1] = 0;
			T_lright[n][1] = 0;
			
			mu_left[n]=0;
			mu_right[n]=0;
			c_pl_left[n]=0;
			c_pl_right[n]=0;
			lamt_left[n]=0;
			lamt_right[n]=0;
			h_fg_lleft[n]=0;
			h_fg_lright[n]=0;
									
			T_v[n][1] = T_c;	//初始温度
			T_sat[n][1] = T_c;
			ro_v[n][1]= ro_sat_given_T(T_v[n][1]);	//气体初始密度为饱和气体的密度
			P_v[n][1]=p_sat_given_T(T_v[n][1]);	//初始压力为饱和压力
			h_fg_v[n]= h_fg_given_T(T_v[n][1]);	//表达式和液体的一样
		}

		if (sort[n]==3)	//左液右气
		{
			v_lleft[n][1] = 0;
			v_lright[n][1] = 0;
			T_lleft[n][1] =T_c;
			T_lright[n][1] = 0;

			mu_left[n]=mu_l_given_T(T_lleft[n][1]);
			mu_right[n]=0;
			c_pl_left[n]=cp_l_given_T(T_lleft[n][1]);
			c_pl_right[n]=0;
			lamt_left[n]=lamt_l_given_T(T_lleft[n][1]);
			lamt_right[n]=0;
			h_fg_lleft[n]= h_fg_given_T(T_lleft[n][1]);
			h_fg_lright[n]=0;
			
			T_v[n][1] = T_c;	//初始温度
			T_sat[n][1] = T_c;
			ro_v[n][1] = ro_sat_given_T(T_v[n][1]);	//气体初始密度为饱和气体的密度
			P_v[n][1] = p_sat_given_T(T_v[n][1]);	//初始压力为饱和压力
			h_fg_v[n] = h_fg_given_T(T_v[n][1]);	//表达式和液体的一样
		}

		if (sort[n]==4)	//左气右液
		{
			v_lleft[n][1]=0;
			v_lright[n][1]=0;
			T_lleft[n][1]=0;
			T_lright[n][1]=T_c;

			mu_left[n]=0;
			mu_right[n]=mu_l_given_T(T_lright[n][1]);
			c_pl_left[n]=0;
			c_pl_right[n]=cp_l_given_T(T_lright[n][1]);
			lamt_left[n]=0;
			lamt_right[n]=lamt_l_given_T(T_lright[n][1]);
			h_fg_lleft[n]=0;
			h_fg_lright[n]= h_fg_given_T(T_lright[n][1]);

			T_v[n][1] = T_c;	//初始温度
			T_sat[n][1] = T_c;
			ro_v[n][1] = ro_sat_given_T(T_v[n][1]);	//气体初始密度为饱和气体的密度
			P_v[n][1] = p_sat_given_T(T_v[n][1]);	//初始压力为饱和压力
			h_fg_v[n] = h_fg_given_T(T_v[n][1]);	//表达式和液体的一样
		}

		T_w[n][1]=T_c;	//壁面温度为冷凝段温度
		P_l[n][2]=0;	//液体初始压力为0Pa
	}
	/* 蒸发段初始温度*/
	for (n=1;n<=length;n=n+1)
	{
		if (heat_sig[n]==1)
		{T_w[n][1]=19.67;}	//蒸发段壁面温度，高于其他部分壁面温度
	}
	/*初始化液体单相对流换热系数*/
	double Re_lleft[length+1],Re_lright[length+1],f,Pr;
	for (n=1;n<=length;n=n+1)
	{
		if (sort[n]==1||sort[n]==3||sort[n]==8)	//1表示液，2表示液-气，8表示液-气-液，计算左侧液体换热系数
		{
			Re_lleft[n]=d*ro_l*fabs(v_lleft[n][1])/mu_left[n];
			Pr = Pr_l_given_T(T_lleft[n][1]);
			if (Re_lleft[n]<=2300)
			{
				h_lleft[n]=4.364*lamt_left[n]/d;
			}
			else if (Re_lleft[n]<=10000)
			{
				f=pow(1.82*log10(Re_lleft[n])-1.64,-2);	//摩擦系数
				h_lleft[n]=lamt_left[n]/d*(f/8)*(Re_lleft[n]-1000)*Pr/(1+12.7*pow(f/8,0.5)*(pow(Pr,2/3)-1));
			}
			else
			{
				if (T_w[n][1]>T_lleft[n][1])
				{
					h_lleft[n]=0.023*lamt_left[n]/d*pow(Re_lleft[n],0.8)*pow(Pr,0.4);
				}
				else
				{
					h_lleft[n]=0.023*lamt_left[n]/d*pow(Re_lleft[n],0.8)*pow(Pr,0.3);
				}
			}
		}
		else
		{
			Re_lleft[n]=0;
			h_lleft[n]=0;
		}
	}
 	for (n=1;n<=length;n=n+1)
	{
		if (sort[n]==1||sort[n]==4||sort[n]==8)	//1表示液，4表示气-液，8表示液-气-液，计算右侧液体单相对流换热系数
		{
			Re_lright[n]=d*ro_l*fabs(v_lright[n][1])/mu_right[n];
			Pr = Pr_l_given_T(T_lright[n][1]);
			if (Re_lright[n]<=2300)
			{
				h_lright[n]=4.364*lamt_right[n]/d;
			}
			else if (Re_lright[n]<=10000)
			{
				f=pow(1.82*log10(Re_lright[n])-1.64,-2);
				h_lright[n]=lamt_right[n]/d*(f/8)*(Re_lright[n]-1000)*Pr/(1+12.7*pow(f/8,0.5)*(pow(Pr,2/3)-1));
			}
			else
			{
				if (T_w[n][1]>T_lright[n][1])
				{
					h_lright[n]=0.023*lamt_right[n]/d*pow(Re_lright[n],0.8)*pow(Pr,0.4);
				}
				else
				{
					h_lright[n]=0.023*lamt_right[n]/d*pow(Re_lright[n],0.8)*pow(Pr,0.3);
				}
			}
		}
		else
		{
			Re_lright[n]=0;
			h_lright[n]=0;
		}
	}

	/*******确定各气泡的状态*********/
	int n_bubble;	//气泡序号，累加变量
	double P_sat;	//气泡的饱和压力，Pa
	double Tr_v[100];	//气泡相变总质量，kg/s，大于0表示液膜蒸发，气泡质量增加
	double V_bubble[100][3];	//气泡体积，m3
	double m_bubble[100][3];	//气泡质量，kg
	double T_vb[100][3];	//气泡温度，K
	double T_sat_vb[100][3];	//气泡饱和温度，K，根据气泡压力算得
	double P_vb[100][3];	//气泡压力,Pa

	if (sort[1]==2||sort[1]==4)	//第1个控制体是气或者气液，意味着气泡跨了坐标原点，接下来要找到该气泡的起点
	{
		n_bubble=1;	//当前气泡序号为1
		while (n_bubble<=n_bubble_total)	//气泡序号小于等于气泡总数
		{
			if (n_bubble==1)	//第1个气泡
			{
				Tr_v[n_bubble]=0;
				n=n_trans[2*n_bubble_total];	//定位到最后一个界面位置，即第1个气泡的起点
				while (n<=length)	//将第1个气泡中坐标原点左边各控制体的相变质量累加
				{
					P_sat=p_sat_given_T(T_w[n][1]);	//壁面温度对应的饱和压力
					T_sat[n][1]=T_sat_given_p(P_v[n][1]);
					if (P_sat>=P_v[n][1])
						{Tr_v_sig[n]=1;}	//蒸发
					else
						{Tr_v_sig[n]=-1;}	//冷凝
					Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][1]*fabs(T_w[n][1]-T_sat[n][1])/h_fg_v[n];
					n=n+1;
				}
				n=1;	//定位到第1个控制体，将第1个气泡中坐标原点右边各控制体的相变质量累加
				while (n<=n_trans[2*n_bubble-1])	//n还没到第1个气泡的终点
				{
					P_sat=p_sat_given_T(T_w[n][1]);
					T_sat[n][1]= T_sat_given_p(P_v[n][1]);
					if (P_sat>=P_v[n][1])
						{Tr_v_sig[n]=1;}	//蒸发
					else
						{Tr_v_sig[n]=-1;}	//冷凝
					Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][1]*fabs(T_w[n][1]-T_sat[n][1])/h_fg_v[n];
					n=n+1;
				}
				/*第1个气泡的体积=          终点控制体体积+               起点控制体体积+                    +中间气控制体体积*/
				V_bubble[n_bubble][1]=A_v*(x_v[n_trans[2*n_bubble-1]][1]+x_v[n_trans[2*n_bubble_total]][1]+(length-1-n_trans[2*n_bubble_total]+n_trans[2*n_bubble-1])*dx);
				m_bubble[n_bubble][1]=ro_v[n_trans[2*n_bubble-1]][1]*V_bubble[n_bubble][1];		//终点控制体的密度*体积		
			}
			else//n_bubble!=1
			{
				Tr_v[n_bubble]=0;
				n=n_trans[2*n_bubble-2];	//定位到气泡起点
				while (n<=n_trans[2*n_bubble-1])	//气泡起点到终点
				{
					P_sat=p_sat_given_T(T_w[n][1]);
					T_sat[n][1]= T_sat_given_p(P_v[n][1]);
					if (P_sat>=P_v[n][1])
						{Tr_v_sig[n]=1;}	//蒸发
					else
						{Tr_v_sig[n]=-1;}	//冷凝
					Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][1]*fabs(T_w[n][1]-T_sat[n][1])/h_fg_v[n];
					n=n+1;
				}
				V_bubble[n_bubble][1]=A_v*(x_v[n_trans[2*n_bubble-2]][1]+x_v[n_trans[2*n_bubble-1]][1]+(n_trans[2*n_bubble-1]-n_trans[2*n_bubble-2]-1)*dx);
				m_bubble[n_bubble][1]=ro_v[n_trans[2*n_bubble-1]][1]*V_bubble[n_bubble][1];
			}
			T_vb[n_bubble][1]=T_v[n_trans[2*n_bubble-1]][1];	//记录第n_bubble个气泡的温度，用终点控制体温度表示
			T_sat_vb[n_bubble][1] = T_sat[n_trans[2 * n_bubble - 1]][1];	//记录第n_bubble个气泡的饱和温度，用终点控制体温度表示
			P_vb[n_bubble][1]=P_v[n_trans[2*n_bubble-1]][1];	//记录第n_bubble个气泡的压力，用终点控制体压力表示
			n_bubble=n_bubble+1;	//依次计算第2，3，，n_bubble_total个气泡的相变质量，体积，质量，温度，压力
		}
	}
	else//第1个控制体不是以气体开始，而是以液体开始，没有气泡跨坐标原点
	{
		n_bubble=1;
		while (n_bubble<=n_bubble_total)
		{
			Tr_v[n_bubble]=0;
			n=n_trans[2*n_bubble-1];	//定位到该气泡的*起点*位置
			while (n<=n_trans[2*n_bubble])	//从起点到终点
				{
					P_sat=p_sat_given_T(T_w[n][1]);
					T_sat[n][1]= T_sat_given_p(P_v[n][1]);
					if (P_sat>=P_v[n][1])
						{Tr_v_sig[n]=1;}	//蒸发
					else
						{Tr_v_sig[n]=-1;}	//冷凝
					Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][1]*fabs(T_w[n][1]-T_sat[n][1])/h_fg_v[n];
					n=n+1;
				}
			V_bubble[n_bubble][1]=A_v*(x_v[n_trans[2*n_bubble-1]][1]+x_v[n_trans[2*n_bubble]][1]+(n_trans[2*n_bubble]-n_trans[2*n_bubble-1]-1)*dx);
			m_bubble[n_bubble][1]=ro_v[n_trans[2*n_bubble-1]][1]*V_bubble[n_bubble][1];	//用起点控制体的密度
			T_vb[n_bubble][1]=T_v[n_trans[2*n_bubble-1]][1];	//用起点控制体的温度表示该气泡的温度
			T_sat_vb[n_bubble][1] = T_sat[n_trans[2 * n_bubble - 1]][1];	//用起点控制体的温度表示该气泡的饱和温度
			P_vb[n_bubble][1]=P_v[n_trans[2*n_bubble-1]][1];	//用起点控制体的压力表示该压力的温度
			n_bubble=n_bubble+1;	//依次计算第2，3，n_bubble_total个气泡的相变质量，体积，质量，温度，压力
		}	
	}

	/******液塞的沸腾蒸发******/
	double Tr_vl_lleft[100], Tr_vl_lright[100];	//液塞左(右)界面沸腾传质速率，kg/s，大于0表示液塞质量减小
	double Tr_vl[100];	//液塞总沸腾传质速率，kg/s，大于0表示液塞质量减小
	double dl_vl_left[100],dl_vl_right[100];	//液塞左(右)界面由于沸腾移动的距离，m，大于0表示液塞长度减小
	int Tr_vl_lleft_sig[length+1],Tr_vl_lright_sig[length+1];	//液体相变标记，1表示沸腾，0表示没有相变
	double Tr_vl_left_each[length+1],Tr_vl_right_each[length+1];	//控制体的左（右）界面沸腾传质速率,kg/s
	int n_liquid;	//液塞序号
	double v_minset=1e-3;	//速度下限10-3 m/s
	/*对每个液塞（最多99个）初始化*/
	for (n=1;n<=99;n=n+1)
	{
		Tr_vl_lleft[n]=0;
		Tr_vl_lright[n]=0;
		Tr_vl[n]=0;
		dl_vl_left[n]=0;
		dl_vl_right[n]=0;
	}
	/*确定各控制体的沸腾标记*/
	for (n=1;n<=length;n=n+1)
	{
		Tr_vl_left_each[n]=0;	//每个控制体左侧液体相变质量速率，kg/s
		Tr_vl_right_each[n]=0;	//每个控制体右侧液体相变质量速率,kg/s
		Tr_vl_lright_sig[n]=0;	//每个控制体左侧液体相变标记初始化
		Tr_vl_lleft_sig[n]=0;	//每个控制体右侧液体相变标记初始化
		h_b[n]=0;	//沸腾换热系数初始化
		if (heat_sig[n]==1)	//表示蒸发段
		{
			if (T_lright[n][1]!=0&&(T_w[n][1]>=T_lright[n][1]))	//壁面温度高于液体温度
			{
				Tr_vl_lright_sig[n]=1;	//沸腾标记
			}
			if (T_lleft[n][1]!=0&&(T_w[n][1]>=T_lleft[n][1]))
			{
				Tr_vl_lleft_sig[n]=1;	//沸腾标记
			}
		}
	}
	
	/****根据液塞左右界面的位置来计算和分配其沸腾质量，计算界面移动距离****/
	/**第1个控制体为气或者气-液**/
	if (sort[1]==2||sort[1]==4)	
	{
		n_liquid=1;	//第1个液塞开始循环
		while (n_liquid<=n_liquid_total)
		{
			/*液塞左界面在蒸发段1*/
			if (n_trans[2*n_liquid-1]>=n_heat_start1&&n_trans[2*n_liquid-1]<=n_heat_end1)
			{
				/*液塞右界面不在蒸发段1*/
				if (n_trans[2*n_liquid]>n_heat_end1)
				{
					n=n_trans[2*n_liquid-1];	//定位到液塞左界面
					/*从左界面到蒸发段1的终点，沸腾相变质量全部累计到左界面*/
					while (n<=n_heat_end1)	
					{
						if (fabs(v_lright[n][1])<=v_minset)	//控制体（右边液体）的速度小于速度下限，沸腾换热系数的速度按下限计算
						{
							h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
						}
						else//控制体速度大于等于下限，则用实际值计算沸腾换热系数
						{
							h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
						}
						if (T_w[n][1]>=T_lright[n][1])//满足沸腾发生的条件
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//液塞在蒸发段1的沸腾相变质量全部累计到左界面
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//每个控制体的沸腾传质量
							if (sort[n]==1)//液控制体的沸腾传质量，左=右
							{
								Tr_vl_left_each[n]=Tr_vl_right_each[n];	
							}
						}
						n=n+1;
					}
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//液塞沸腾总质量速率
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//液塞左界面的移动距离
					/*液塞右界面在蒸发段2，从蒸发段2的起点到液塞右界面沸腾传质量累加到右界面*/
					if ((n_trans[2*n_liquid]>=n_heat_start2)&&(n_trans[2*n_liquid]<=n_heat_end2))
					{
						n=n_heat_start2;	//定位到蒸发段2的起点
						while (n<=n_trans[2*n_liquid]-1)
						{
							if (fabs(v_lright[n][1])<=v_minset)//用速度计算沸腾换热系数
							{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
							else
							{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
							if (T_w[n][1]>=T_lright[n][1])//满足沸腾条件
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质量累加到右界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//每个控制体右界面沸腾传质量
								if (sort[n]==1)
							    {Tr_vl_left_each[n]=Tr_vl_right_each[n];}	//液控制体沸腾传质量，左=右
							}
							n=n+1;
						}
						n=n_trans[2*n_liquid];	//定位到液塞右界面所在控制体
						if (fabs(v_lleft[n][1])<=v_minset)
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
						else
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//与前面不同之处在于x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//与前面不同之处在于x_lleft[n][1]
							if (sort[n]==1)
							{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
						}
						if (n_liquid==n_liquid_total)//液塞右界面的沸腾总质量累加到下一个液塞？？？
						{Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];}
						else
						{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];}
					
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面由于沸腾产生的移动距离
					}
					/*液塞右界面超过了蒸发段2，蒸发段2的沸腾质量累加到右界面*/
					if (n_trans[2*n_liquid]>n_heat_end2)
					{
						n=n_heat_start2;
						while (n<=n_heat_end2)
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
							else
							{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
							}
							n=n+1;
						}
						if (n_liquid==n_liquid_total)//液塞右界面的沸腾总质量累加到下一个液塞？？？
						{Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];}
						else
						{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];}
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);
					}
				}
				/*液塞右界面在蒸发段1，从左界面到右界面的所有控制体，沸腾质量累加到左界面，然后平均分给两个相邻液塞*/
				else
				{
					n=n_trans[2*n_liquid-1];	//定位到液塞左界面
					while (n<=(n_trans[2*n_liquid]-1))	
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;
					}
					n=n_trans[2*n_liquid];	//定位到右界面所在控制体
					if (fabs(v_lleft[n][1])<=v_minset)
					{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
					else
					{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
					if (T_w[n][1]>=T_lleft[n][1])
					{
						Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft
						Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
						if (sort[n]==1)
						{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
					}
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//分给n_liquid一半
					if (n_liquid==n_liquid_total)	//分给n_liquid+1一半
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lleft[n_liquid]/2;}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lleft[n_liquid]/2;}
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//只有左界面位移？？？
				}	
			}
			/*液塞左界面在蒸发段2*/
			else if (n_trans[2*n_liquid-1]>=n_heat_start2&&n_trans[2*n_liquid-1]<=n_heat_end2)
			{
				/*液塞右界面超过了蒸发段2*/
				if (n_trans[2*n_liquid]>n_heat_end2)
				{
					n=n_trans[2*n_liquid-1];	//定位到左界面
					/*从左界面到蒸发段2的终点，沸腾相变质量全部累计到左界面*/
					while (n<=n_heat_end2)
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;
					}
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//液塞沸腾总质量
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//液塞左界面的移动距离
				}
				/*液塞右界面没超过蒸发段2，则从左界面到右界面的所有控制体，沸腾质量累加到左界面，然后平均分给两个相邻液塞？？？*/
				else
				{
					n=n_trans[2*n_liquid-1];	//定位到左界面
					while (n<=(n_trans[2*n_liquid]-1))	//从左界面到右界面左边的控制体，沸腾质量累加
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;
					}
					n=n_trans[2*n_liquid];	//定位到右界面
					if (fabs(v_lleft[n][1])<=v_minset)
					{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
					else
					{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
					if (T_w[n][1]>=T_lleft[n][1])
					{
						Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
						Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
						if (sort[n]==1)
						{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
					}
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//n_liquid分配一半
					if (n_liquid==n_liquid_total)	//n_liquid+1分配一半
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lleft[n_liquid]/2;}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lleft[n_liquid]/2;}
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//位移全部分配给左界面
				}
			}
			/*左界面在1到蒸发段1的起点之间*/
			else if (n_trans[2*n_liquid-1]>=1&&n_trans[2*n_liquid-1]<n_heat_start1)
			{
				/*右界面在蒸发段1*/
				if (n_trans[2*n_liquid]>=n_heat_start1&&n_trans[2*n_liquid]<=n_heat_end1)
				{
					n=n_heat_start1;	//定位到蒸发段1起点
					/*从蒸发段1起点到液塞右界面，沸腾质量累加到右界面*/
					while (n<=n_trans[2*n_liquid]-1)
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;	
					}
					n=n_trans[2*n_liquid];	//定位到右界面
					if (fabs(v_lleft[n][1])<=v_minset)
					{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
					else
					{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
					if (T_w[n][1]>=T_lleft[n][1])
					{
						Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
						Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
						if (sort[n]==1)
						{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
					}
					/*右界面的沸腾质量分配给下一个液塞？？*/
					if (n_liquid==n_liquid_total)
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];}
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);
				}
				/*右界面在蒸发段1之后，在蒸发段2之前*/
				else if (n_trans[2*n_liquid]>n_heat_end1&&n_trans[2*n_liquid]<n_heat_start2)
				{
					n=n_heat_start1;	//定位到蒸发段1起点
					/*蒸发段1的沸腾质量累加到左界面，然后平均分给两个相邻液塞*/
					while (n<=n_heat_end1)
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;
					}
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//n_liquid分配一半
					if (n_liquid==n_liquid_total)	//n_liquid+1分配一半
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lleft[n_liquid]/2;}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lleft[n_liquid]/2;}
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面的移动距离
				}
				/*右界面在蒸发段2*/
				else if (n_trans[2*n_liquid]>=n_heat_start2&&n_trans[2*n_liquid]<=n_heat_end2)
				{
					n=n_heat_start1;	//定位到蒸发段1起点
					/*蒸发段1的沸腾质量累加到左界面*/
					while (n<=n_heat_end1)
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;
					}
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//沸腾质量分配给n_liquid
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
					n=n_heat_start2;	//定位到蒸发段2起点
					/*蒸发段2起点到右界面的沸腾质量累加到右界面*/
					while (n<=n_trans[2*n_liquid])
					{
						if (fabs(v_lleft[n][1])<=v_minset)
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
						else
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
						}
						n=n+1;
					}
					if (n_liquid==n_liquid_total)	//沸腾质量分配给下一个液塞
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];}
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
				}
				/*右界面超过蒸发段2*/
				else if (n_trans[2*n_liquid]>n_heat_end2)
				{
					n=n_heat_start1;	//定位到蒸发段1起点
					/*蒸发段1沸腾质量累加到左界面*/
					while (n<=n_heat_end1)
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;
					}
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//沸腾质量分配给n_liquid
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
					n=n_heat_start2;	//定位到蒸发段2起点
					/*蒸发段2沸腾质量累加到右界面*/
					while (n<=n_heat_end2)
					{
						if (fabs(v_lleft[n][1])<=v_minset)
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
						else
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
						}
						n=n+1;
					}
					if (n_liquid==n_liquid_total)	//沸腾质量分配给n_liquid+1
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];}
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
				}

			}
			/*左界面在蒸发段1之后，在蒸发段2之前*/
			else if (n_trans[2*n_liquid-1]>n_heat_end1&&n_trans[2*n_liquid-1]<n_heat_start2)
			{
				/*右界面在蒸发段2*/
				if (n_trans[2*n_liquid]>=n_heat_start2&&n_trans[2*n_liquid]<=n_heat_end2)
				{
					n=n_heat_start2;	//定位到蒸发段2起点
					/*蒸发段2起点到右界面的沸腾质量累加到右界面*/
					while (n<=n_trans[2*n_liquid])
					{
						if (fabs(v_lleft[n][1])<=v_minset)
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
						else
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
						}
						n=n+1;
					}
					if (n_liquid==n_liquid_total)	//沸腾质量分配到下一个液塞
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];}
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
				}
				/*右界面超过蒸发段2*/
				else if (n_trans[2*n_liquid]>n_heat_end2)
				{
					n=n_heat_start2;	//定位到蒸发段2起点
					/*蒸发段2的沸腾质量累加到右界面，平均分配给相邻两个液塞*/
					while (n<=n_heat_end2)
					{
						if (fabs(v_lleft[n][1])<=v_minset)
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
						else
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
						}
						n=n+1;
					}
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid]/2;	//n_liquid分配一半
					if (n_liquid==n_liquid_total)	//n_liquid+1分配一半
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid]/2;}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid]/2;}
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面的移动距离
				}
			}
			/*液塞左右界面位置分类讨论完毕，考虑下一个液塞*/
			n_liquid=n_liquid+1;
		}
	}
	/**第1个控制体为液sort[1]==1或者液-气sort[1]==3或者液-气-液sort[1]==8**/
	else                                
	{
		n_liquid=1;   //第1个液塞单独考虑，是因为坐标原点将第1个液塞分成了两段，先考虑原点右边这一段，再考虑原点左边那一段
		/*第1个液塞的右界面（n_trans[1]）在蒸发段1*/
		if (n_trans[1]>=n_heat_start1&&n_trans[1]<=n_heat_end1)
		{
			n=n_heat_start1;	//定位到蒸发段1起点
			/*蒸发段1起点到右界面的沸腾质量累加到右界面*/
			while (n<=n_trans[1])
			{
				if (fabs(v_lleft[n][1])<=v_minset)
				{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
				else
				{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
				if (T_w[n][1]>=T_lleft[n][1])
				{
					Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
					Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
					if (sort[n]==1)
					{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
				}
				n=n+1;
			}
			Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//沸腾质量分配给n_liquid
			dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
		}
		/*第1个液塞右界面在蒸发段1之后，蒸发段2之前*/
		else if (n_trans[1]>n_heat_end1&&n_trans[1]<n_heat_start2)
		{
			n=n_heat_start1;	//定位到蒸发段1起点
			/*蒸发段1沸腾质量累加到右界面，平均分配给相邻两液塞*/
			while (n<=n_heat_end1)
			{
				if (fabs(v_lleft[n][1])<=v_minset)
				{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
				else
				{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
				if (T_w[n][1]>=T_lleft[n][1])
				{
					Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
					Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
					if (sort[n]==1)
					{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
				}
				n=n+1;
			}
			Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid]/2;	//n_liquid分配一半
			Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lright[n_liquid]/2;	//n_liquid_total分配一半
			dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
		}
		/*第1个液塞右界面在蒸发段2*/
		else if (n_trans[1]>=n_heat_start2&&n_trans[1]<=n_heat_end2)
		{
			n=n_heat_start1;	//定位到蒸发段1起点
			/*蒸发的1的沸腾质量累加到左界面？？？*/
			while (n<=n_heat_end1)
			{
				if (fabs(v_lleft[n][1])<=v_minset)
				{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
				else
				{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
				if (T_w[n][1]>=T_lleft[n][1])
				{
					Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
					Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
					if (sort[n]==1)
					{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
				}
				n=n+1;
			}
			Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid];	//沸腾质量分配给n_liquid_total
			dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
			n=n_heat_start2;	//定位到蒸发段2起点
			/*蒸发段2起点到右界面的沸腾质量累加到右界面*/
			while (n<=n_trans[1])
			{
				if (fabs(v_lleft[n][1])<=v_minset)
				{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
				else
				{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
				if (T_w[n][1]>=T_lleft[n][1])
				{
					Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
					Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
					if (sort[n]==1)
					{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
				}
				n=n+1;
			}
			Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//沸腾质量分配给n_liquid
			dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
		}

		/*第1个液塞左界面(n_trans[2*n_liquid_total]）在蒸发段1*/
		if (n_trans[2*n_liquid_total]>=n_heat_start1&&n_trans[2*n_liquid_total]<=n_heat_end1)
		{
			n=n_trans[2*n_liquid_total];	//定位到左界面
			/*左界面到蒸发段1终点沸腾质量累加到左界面*/
			while (n<=n_heat_end1)
			{
				if (fabs(v_lright[n][1])<=v_minset)
				{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
				else
				{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
				if (T_w[n][1]>=T_lright[n][1])
				{
					Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
					Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
					if (sort[n]==1)
					{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
				}
				n=n+1;
			}
			Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid];	//沸腾质量分配给n_liquid_total
			dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
			n=n_heat_start2;	//定位到蒸发段2起点
			/*蒸发段2沸腾质量累加到右界面*/
			while (n<=n_heat_end2)
			{
				if (fabs(v_lright[n][1])<=v_minset)
				{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
				else
				{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
				if (T_w[n][1]>=T_lright[n][1])
				{
					Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];		//x_lright[n][1]
					Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
					if (sort[n]==1)
					{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
				}
				n=n+1;
			}
			Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//沸腾质量分配给n_liquid=1
			dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
		}
		/*第1个液塞左界面在蒸发段1之后，在蒸发段2之前*/
		else if (n_trans[2*n_liquid_total]>n_heat_end1&&n_trans[2*n_liquid_total]<n_heat_start2)
		{
			n=n_heat_start2;	//定位到蒸发段2起点
			/*蒸发段2的沸腾质量累加到左界面，平均分配给相邻两液塞*/
			while (n<=n_heat_end2)
			{
				if (fabs(v_lright[n][1])<=v_minset)
				{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
				else
				{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
				if (T_w[n][1]>=T_lright[n][1])
				{
					Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
					Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
					if (sort[n]==1)
					{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
				}
				n=n+1;
			}
			Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid]/2;	//分配给n_liquid_total一半
			Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//分配给n_liquid一半
			dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
		}
		/*第1个液塞左界面在蒸发段2*/
		else if (n_trans[2*n_liquid_total]>=n_heat_start2&&n_trans[2*n_liquid_total]<=n_heat_end2)
		{
			n=n_trans[2*n_liquid_total];	//定位到左界面
			/*左界面到蒸发段2终点的沸腾质量累加到左界面*/
			while (n<=n_heat_end2)
			{
				if (fabs(v_lright[n][1])<=v_minset)
				{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
				else
				{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
				if (T_w[n][1]>=T_lright[n][1])
				{
					Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
					Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
					if (sort[n]==1)
					{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
				}
				n=n+1;
			}
			Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid];	//沸腾质量分配给n_liquid_total
			dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
		}

		n_liquid=2;	//第2个液塞开始循环
		while (n_liquid<=n_liquid_total)
		{
			/*左界面n_trans[2*n_liquid-2]在蒸发段1*/
			if (n_trans[2*n_liquid-2]>=n_heat_start1&&n_trans[2*n_liquid-2]<=n_heat_end1)
			{
				/*右界面超过蒸发段1*/
				if (n_trans[2*n_liquid-1]>n_heat_end1)
				{
					n=n_trans[2*n_liquid-2];	//定位到左界面
					/*左界面到蒸发段1终点的沸腾质量累加到左界面*/
					while (n<=n_heat_end1)
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;
					}
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//沸腾质量分配给n_liquid-1
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
					/*右界面在蒸发段2*/
					if ((n_trans[2*n_liquid-1]>=n_heat_start2)&&(n_trans[2*n_liquid-1]<=n_heat_end2))
					{
						n=n_heat_start2;	//定位到蒸发段2起点
						/*蒸发段2起点到右界面的沸腾质量累加到右界面*/
						while (n<=n_trans[2*n_liquid-1]-1)
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
							else
							{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
							}
							n=n+1;
						}
						n=n_trans[2*n_liquid-1];	//定位到右界面控制体，单独计算沸腾质量并累加到右界面
						if (fabs(v_lleft[n][1])<=v_minset)
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
						else
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//沸腾质量分配给n_liquid
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
					}
					/*右界面超过蒸发段2*/
					if (n_trans[2*n_liquid-1]>n_heat_end2)
					{
						n=n_heat_start2;	//定位到蒸发段2起点
						/*蒸发段2沸腾质量累加到右界面*/
						while (n<=n_heat_end2)
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
							else
							{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//沸腾质量分配给n_liquid
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
					}
				}
				/*右界面在蒸发段1*/
				else
				{
					n=n_trans[2*n_liquid-2];	//定位到左界面
					/*从左界面到右界面沸腾质量累加到左界面，平均分配到两个相邻气泡*/
					while (n<=(n_trans[2*n_liquid-1]-1))
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;
					}
					n=n_trans[2*n_liquid-1];	//定位到右界面，单独计算沸腾质量并累加到左界面
					if (fabs(v_lleft[n][1])<=v_minset)
					{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
					else
					{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
					if (T_w[n][1]>=T_lleft[n][1])
					{
						Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
						Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
						if (sort[n]==1)
						{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
					}
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid]/2;	//分配一半给n_liquid-1
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//分配一半给n_liquid
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
				}	
			}
			/*左界面在蒸发段2*/
			else if (n_trans[2*n_liquid-2]>=n_heat_start2&&n_trans[2*n_liquid-2]<=n_heat_end2)
			{
				/*右界面超过蒸发段2*/
				if (n_trans[2*n_liquid-1]>n_heat_end2)
				{
					n=n_trans[2*n_liquid-2];	//定位到左界面
					/*左界面到蒸发段2终点的沸腾质量累加到左界面*/
					while (n<=n_heat_end2)
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;
					}
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//沸腾质量分配给n_liquid-1
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
				}
				/*右界面在蒸发段2*/
				else
				{
					n=n_trans[2*n_liquid-2];	//定位到左界面
					/*从左界面到右界面沸腾质量累加到左界面，平均分配到两个相邻气泡*/
					while (n<=(n_trans[2*n_liquid-1]-1))
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;
					}
					n=n_trans[2*n_liquid-1];	//定位到右界面控制体，单独计算沸腾质量，累加到左界面
					if (fabs(v_lleft[n][1])<=v_minset)
					{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
					else
					{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
					if (T_w[n][1]>=T_lleft[n][1])
					{
						Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
						Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
						if (sort[n]==1)
						{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
					}
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid]/2;	//分配一半给n_liquid-1
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//分配一半给n_liquid
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
				}
			}
			/*左界面在1到蒸发段1起点之间*/
			else if (n_trans[2*n_liquid-2]>=1&&n_trans[2*n_liquid-2]<n_heat_start1)
			{
				/*右界面在蒸发段1*/
				if (n_trans[2*n_liquid-1]>=n_heat_start1&&n_trans[2*n_liquid-1]<=n_heat_end1)
				{
					n=n_heat_start1;	//定位到蒸发段1起点
					/*蒸发段1起点到右界面沸腾质量累加到右界面*/
					while (n<=n_trans[2*n_liquid-1]-1)
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;	
					}
					n=n_trans[2*n_liquid-1];	//定位到右界面控制体，单独计算沸腾质量，并累加到右界面
					if (fabs(v_lleft[n][1])<=v_minset)
					{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
					else
					{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
					if (T_w[n][1]>=T_lleft[n][1])
					{
						Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
						Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
						if (sort[n]==1)
						{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
					}
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//分配给n_liquid
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
				}
				/*右界面在蒸发段1之后，蒸发段2之前*/
				else if (n_trans[2*n_liquid-1]>n_heat_end1&&n_trans[2*n_liquid-1]<n_heat_start2)
				{
					n=n_heat_start1;	//定位到蒸发段1起点
					/*蒸发段1的沸腾质量累加到左界面，平均分配给两个相邻气泡*/
					while (n<=n_heat_end1)
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;
					}
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid]/2;	//分配一半给n_liquid-1，即左侧气泡
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//分配一半给n_liquid，即右侧气泡
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左侧界面移动距离
				}
				/*右界面在蒸发段2*/
				else if (n_trans[2*n_liquid-1]>=n_heat_start2&&n_trans[2*n_liquid-1]<=n_heat_end2)
				{
					n=n_heat_start1;	//定位到蒸发段1起点
					/*蒸发段1沸腾质量累加到左界面*/
					while (n<=n_heat_end1)
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;
					}
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//分配给n_liquid-1，即左侧气泡
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
					n=n_heat_start2;	//定位到蒸发段2起点
					/*蒸发段2起点到右界面沸腾质量累加到右界面*/
					while (n<=n_trans[2*n_liquid-1])
					{
						if (fabs(v_lleft[n][1])<=v_minset)
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
						else
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
						}
						n=n+1;
					}
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//沸腾质量分配给n_liquid，即右侧气泡
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
				}
				/*右界面超过蒸发段2*/
				else if (n_trans[2*n_liquid-1]>n_heat_end2)
				{
					n=n_heat_start1;	//定位到蒸发段1起点
					/*蒸发段1沸腾质量累加到左界面*/
					while (n<=n_heat_end1)
					{
						if (fabs(v_lright[n][1])<=v_minset)
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
						else
						{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
						if (T_w[n][1]>=T_lright[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//x_lright[n][1]
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
							if (sort[n]==1)
							{Tr_vl_left_each[n]=Tr_vl_right_each[n];}
						}
						n=n+1;
					}
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//沸腾质量分配给n_liquid-1，即左侧气泡
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
					n=n_heat_start2;	//定位到蒸发段2起点
					/*蒸发段2沸腾质量累加到右界面*/
					while (n<=n_heat_end2)
					{
						if (fabs(v_lleft[n][1])<=v_minset)
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
						else
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
						}
						n=n+1;
					}
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//沸腾质量分配给n_liquid，即右侧气泡
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
				}
			}
			/*左界面在蒸发段1之后，蒸发段2之前*/
			else if (n_trans[2*n_liquid-2]>n_heat_end1&&n_trans[2*n_liquid-2]<n_heat_start2)
			{
				/*右界面在蒸发段2*/
				if (n_trans[2*n_liquid-1]>=n_heat_start2&&n_trans[2*n_liquid-1]<=n_heat_end2)
				{
					n=n_heat_start2;	//定位到蒸发段2起点
					/*蒸发段2沸腾质量累加到右界面*/
					while (n<=n_trans[2*n_liquid-1])
					{
						if (fabs(v_lleft[n][1])<=v_minset)
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
						else
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
						}
						n=n+1;
					}
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//沸腾质量分配给n_liquid，即右侧气泡
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
				}
				/*右界面超过蒸发段2*/
				else if (n_trans[2*n_liquid-1]>n_heat_end2)
				{
					n=n_heat_start2;	//定位到蒸发段2起点
					/*蒸发段2沸腾质量累加到右界面，平均分配给两个相邻气泡*/
					while (n<=n_heat_end2)
					{
						if (fabs(v_lleft[n][1])<=v_minset)
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
						else
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
						}
						n=n+1;
					}
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lright[n_liquid]/2;	//分配右边给n_liquid-1，即左侧气泡
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid]/2;	//分配一半给n_liquid，即右侧气泡
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
				}
			}
			/*分类讨论完毕，下一个液塞*/
			n_liquid=n_liquid+1;
		}

	}

	printf("初始化完成\n");	//屏幕输出

	FILE *fid0;	//定义一个叫fid0的文件指针
	fid0 =fopen("state.txt","a+");	//按照a+模式打开一个名叫state.txt的文件，返回state.txt的文件指针给fid0,异常则返回NULL
	/*a+表示以附加方式打开可读写的文件。若文件不存在，则会建立该文件，如果文件存在，写入的数据会被加到文件尾后，
	即文件原先的内容会被保留。原来的EOF符不保留*/
	fprintf(fid0,"初始化完成\n");	//将“初始化完成\n”写入由fid0指出的文件，格式与printf类似，可以指明写入数据类型
	fclose(fid0);	//关闭由fid0指出的文件,返回操作结果，0或EOF

	FILE *fid1=fopen("ro_v.txt","a+");	//按照a+模式打开一个名叫ro_v.txt的文件，返回ro_v.txt的文件指针给fid1
	FILE *fid2=fopen("T_w.txt","a+");	//按照a+模式打开一个名叫T_w.txt的文件，返回T_w.txt的文件指针给fid2
	FILE *fid3=fopen("trans_point.txt","a+");
	FILE *fid4=fopen("P_v.txt","a+");
	FILE *fid5=fopen("T_e.txt","a+");
	FILE *fid6=fopen("x_l.txt","a+");
	FILE *fid7=fopen("T_e2.txt","a+");
	FILE *fid8=fopen("T_v.txt","a+");
	FILE *fid9=fopen("v_l.txt","a+");
	FILE *fid10=fopen("T_l.txt","a+");
	FILE *fid11=fopen("sensible.txt","a+");
	FILE *fid12=fopen("latent.txt","a+");
	FILE *fid13=fopen("T_vb.txt","a+");
	FILE *fid14=fopen("P_vb.txt","a+");
	FILE *fid15=fopen("ro_vb.txt","a+");
	FILE *fid16=fopen("P_l.txt","a+");
	FILE *fid17=fopen("liquid_total.txt","a+");
	FILE *fid18=fopen("bubble_gen_total.txt","a+");
	FILE *fid19=fopen("P_place.txt","a+");
	FILE *fid20=fopen("h_b.txt","a+");
	FILE *fid21=fopen("h_l.txt","a+");
	FILE *fid22=fopen("T_sat.txt", "a+");
	FILE *fid23 = fopen("T_sat_vb.txt", "a+");
	
	/*写入第一行表头表示n从1到length，然后换行*/
	int location = 1;       //输出时的数据坐标，相当于n*/
	for (location = 1; location <= length; location++)
	{
		fprintf(fid1, " %d", location);	//空格用于分列，n=1在第2列
		fprintf(fid2, " %d", location);	//空格用于分列，n=1在第2列
		fprintf(fid4, " %d", location);	//空格用于分列，n=1在第2列
		fprintf(fid6, " %d", location);	//空格用于分列，n=1在第2列
		fprintf(fid8, " %d", location);	//空格用于分列，n=1在第2列
		fprintf(fid9, " %d", location);	//空格用于分列，n=1在第2列
		fprintf(fid10, " %d", location);	//空格用于分列，n=1在第2列
		fprintf(fid16, " %d", location);	//空格用于分列，n=1在第2列
		fprintf(fid20, " %d", location);	//空格用于分列，n=1在第2列
		fprintf(fid21, " %d", location);	//空格用于分列，n=1在第2列
		fprintf(fid22, " %d", location);	//空格用于分列，n=1在第2列
	}
	/*表头遍历n之后换行*/
	fprintf(fid1, "\n");
	fprintf(fid2, "\n");
	fprintf(fid4, "\n");
	fprintf(fid6, "\n");
	fprintf(fid8, "\n");
	fprintf(fid9, "\n");
	fprintf(fid10, "\n");
	fprintf(fid16, "\n");
	fprintf(fid20, "\n");
	fprintf(fid21, "\n");
	fprintf(fid22, "\n");

	fclose(fid1);	//关闭由fid1指出的文件
	fclose(fid2);	//关闭由fid2指出的文件
	fclose(fid3);
	fclose(fid4);
	fclose(fid5);
	fclose(fid6);
	fclose(fid7);
	fclose(fid8);
	fclose(fid9);
	fclose(fid10);
	fclose(fid11);
	fclose(fid12);
	fclose(fid13);
	fclose(fid14);
	fclose(fid15);
	fclose(fid16);
	fclose(fid17);
	fclose(fid18);
	fclose(fid19);
	fclose(fid20);
	fclose(fid21);
	fclose(fid22);
	fclose(fid23);

	/*初始化下一时间节点的值*/
	for (n=1;n<=length;n=n+1)
	{
		x_lleft[n][2]=0;	//控制体左边液塞长度
		x_lright[n][2]=0;	//控制体右边液塞长度
		v_lleft[n][2]=0;	//控制体左边液塞速度
		v_lright[n][2]=0;
		T_lleft[n][2]=0;	//控制体左边液塞温度
		T_lright[n][2]=0;
		x_v[n][2]=0;	//控制体中气塞长度
		ro_v[n][2]=0;	//气塞密度
		T_v[n][2]=0;	//气塞温度
		T_sat[n][2] = 0;	//气塞饱和温度
		P_v[n][2]=0;	//气塞压力
	}

	/*******开始主循环*********/

    int i=1;	//i是时间节点循环控制变量
	int kk=0;	//输出的数据组数
	double sensible_c_total=0;	//冷凝段总显热传热量（多个时间步长下的累加量），W
	double sensible_e_total=0;	//蒸发段总显热传热量（多个时间步长下的累加量），W
	double latent_c_total=0;	//冷凝段总潜热传热量（多个时间步长下的累加量），W
	double latent_e_total=0;	//蒸发段总潜热传热量（多个时间步长下的累加量），W

	while (i<=nn)
	{
		int cross_positive=0;	//液塞左界面移动标记，=1表示正向跨过原点，第1个气泡变第2个 
		int cross_negative=0;	//液塞左界面移动标记，=1表示负向跨过原点，第1个气泡变最后一个
		double sensible_c=0;	//冷凝段显热（液体单相对流换热），W
		double latent_c=0;	//冷凝段潜热,W
		double sensible_e=0;	//蒸发段显热,W
		double latent_e=0;	//蒸发段潜热,W

		/*******计算传热量（显热和潜热）*********/
		n=1;
		while (n<=length)
		{
			if (heat_sig[n]==-1)	//-1表示冷凝
			{
				/*计算冷凝段显热和潜热*/
				if (sort[n]==1)	//1表示液控制体
				{
					sensible_c=sensible_c+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d;	//单相对流，控制体左边的液体温度长度换热系数来表示
				}
				else //其他含有气体的控制体，2气，3液-气，4气-液，8液-气-液
				{
					sensible_c=sensible_c+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d;	//单相对流，左边液体温度长度换热系数
					sensible_c=sensible_c+h_lright[n]*(T_w[n][1]-T_lright[n][1])*x_lright[n][1]*pi*d;	//单相对流，右边液体温度长度换热系数
					latent_c=latent_c+Tr_v_sig[n]*h_v[n]*fabs(T_w[n][1]-T_sat[n][1])*x_v[n][1]*pi*d;	//液膜的蒸发/冷凝换热
				}
			}
			if (heat_sig[n]==1)	//1表示蒸发段
			{
				/*计算蒸发段显热和潜热*/
				if (sort[n]==1)		//液控制体			
				{
					sensible_e=sensible_e+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d;	//单相对流，左边液体参数
					latent_e=latent_e+Tr_vl_lleft_sig[n]*h_b[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d;	//沸腾，左边液体参数
				}
				else//含气控制体
				{
					sensible_e=sensible_e+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d;	//单相对流，左
					sensible_e=sensible_e+h_lright[n]*(T_w[n][1]-T_lright[n][1])*x_lright[n][1]*pi*d;	//单相对流，右
					latent_e=latent_e+Tr_v_sig[n]*h_v[n]*fabs(T_w[n][1]-T_sat[n][1])*x_v[n][1]*pi*d
							+Tr_vl_lleft_sig[n]*h_b[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
							+Tr_vl_lright_sig[n]*h_b[n]*(T_w[n][1]-T_lright[n][1])*x_lright[n][1]*pi*d;	//液膜蒸发+左沸腾+右沸腾
				}
			}
			n=n+1;
		}

		/*******计算壁面温度（第1个和最后1个控制体单独考虑）*******/
		n=1;
		if (plain[n]==0)	//0表示壁温可变
		{
			if (sort[n]==1)	//液控制体
			{
				T_w[n][2]=T_w[n][1]+dt/(c_w*ro_w*A_w*dx)*(heat_flux[n]*dx-h_c[n]*pi*d_out*(T_w[n][1]-T_out)*dx
						-Tr_v_sig[n]*h_v[n]*fabs(T_w[n][1]-T_sat[n][1])*x_v[n][1]*pi*d  //这一项多余？？？
						-Tr_vl_lleft_sig[n]*h_b[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
						-h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
						-lamt_w*A_w*(T_w[n][1]-T_w[length][1])/dx-lamt_w*A_w*(T_w[n][1]-T_w[n+1][1])/dx);	//壁面能量方程
			}
			else//含气控制体
			{
				T_w[n][2]=T_w[n][1]+dt/(c_w*ro_w*A_w*dx)*(heat_flux[n]*dx-h_c[n]*pi*d_out*(T_w[n][1]-T_out)*dx
						-Tr_v_sig[n]*h_v[n]*fabs(T_w[n][1]-T_sat[n][1])*x_v[n][1]*pi*d
						-Tr_vl_lleft_sig[n]*h_b[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
						-Tr_vl_lright_sig[n]*h_b[n]*(T_w[n][1]-T_lright[n][1])*x_lright[n][1]*pi*d
						-h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
						-h_lright[n]*(T_w[n][1]-T_lright[n][1])*x_lright[n][1]*pi*d
						-lamt_w*A_w*(T_w[n][1]-T_w[length][1])/dx-lamt_w*A_w*(T_w[n][1]-T_w[n+1][1])/dx);
			}
		}
		else//壁温已被设定，则为初始值，只有冷凝段
			{T_w[n][2]=T_w[n][1];}

		n=2;
		while (n<=(length-1))
		{
			if (plain[n]==0)
			{
				if (sort[n]==1)
				{
					T_w[n][2]=T_w[n][1]+dt/(c_w*ro_w*A_w*dx)*(heat_flux[n]*dx-h_c[n]*pi*d_out*(T_w[n][1]-T_out)*dx
						-Tr_v_sig[n]*h_v[n]*fabs(T_w[n][1]-T_sat[n][1])*x_v[n][1]*pi*d
						-Tr_vl_lleft_sig[n]*h_b[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
						-h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
						-lamt_w*A_w*(T_w[n][1]-T_w[n-1][1])/dx-lamt_w*A_w*(T_w[n][1]-T_w[n+1][1])/dx);
				}
				else
				{
					T_w[n][2]=T_w[n][1]+dt/(c_w*ro_w*A_w*dx)*(heat_flux[n]*dx-h_c[n]*pi*d_out*(T_w[n][1]-T_out)*dx
						-Tr_v_sig[n]*h_v[n]*fabs(T_w[n][1]-T_sat[n][1])*x_v[n][1]*pi*d
						-Tr_vl_lleft_sig[n]*h_b[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
						-Tr_vl_lright_sig[n]*h_b[n]*(T_w[n][1]-T_lright[n][1])*x_lright[n][1]*pi*d
						-h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
						-h_lright[n]*(T_w[n][1]-T_lright[n][1])*x_lright[n][1]*pi*d
						-lamt_w*A_w*(T_w[n][1]-T_w[n-1][1])/dx-lamt_w*A_w*(T_w[n][1]-T_w[n+1][1])/dx);
				}
			}
		else
			{T_w[n][2]=T_w[n][1];}
			n=n+1;
		}

		n=length;
		if (plain[n]==0)
		{
			if (sort[n]==1)
			{
				T_w[n][2]=T_w[n][1]+dt/(c_w*ro_w*A_w*dx)*(heat_flux[n]*dx-h_c[n]*pi*d_out*(T_w[n][1]-T_out)*dx
					-Tr_v_sig[n]*h_v[n]*fabs(T_w[n][1]-T_sat[n][1])*x_v[n][1]*pi*d
					-Tr_vl_lleft_sig[n]*h_b[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
					-h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
					-lamt_w*A_w*(T_w[n][1]-T_w[n-1][1])/dx-lamt_w*A_w*(T_w[n][1]-T_w[1][1])/dx);
			}
			else
			{
				T_w[n][2]=T_w[n][1]+dt/(c_w*ro_w*A_w*dx)*(heat_flux[n]*dx-h_c[n]*pi*d_out*(T_w[n][1]-T_out)*dx
					-Tr_v_sig[n]*h_v[n]*fabs(T_w[n][1]-T_sat[n][1])*x_v[n][1]*pi*d
					-Tr_vl_lleft_sig[n]*h_b[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
					-Tr_vl_lright_sig[n]*h_b[n]*(T_w[n][1]-T_lright[n][1])*x_lright[n][1]*pi*d
					-h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
					-h_lright[n]*(T_w[n][1]-T_lright[n][1])*x_lright[n][1]*pi*d
					-lamt_w*A_w*(T_w[n][1]-T_w[n-1][1])/dx-lamt_w*A_w*(T_w[n][1]-T_w[1][1])/dx);
			}
		}
		else
			{T_w[n][2]=T_w[n][1];}

		/*******初始化每个控制体的相变标记*******/
		for (n=1;n<=length;n=n+1)
			{Tr_v_sig[n]=0;}//0表示无相变

		/*******定义一些量*******/
		int jump;	//左界面移动标记，0表示没有移出原控制体，1表示向右移出原控制体
		double dm_left[100], dm_right[100];	//相邻气泡相变导致的液塞左（右）界面净流入质量,大于0表示液塞质量增加
		double dl_left[100], dl_right[100];	//相邻气泡相变导致的液塞左（右)界面总移动距离，大于0表示液塞长度增加
		double G_force[100];	//液塞重力，N/m2
		double l[100];	//液塞长度？
		double C_ll;	//液塞总阻力系数，大于0表示向右，小于0表示向左
		double K;	//弯头阻力系数

		/*******对于起始点（n=1）状态的分类讨论：第一部分：气或者气-液，此时液塞序号领先于气泡序号*********/
		/*******计算各液塞的质量变化和重力作用、液塞总长度、速度、温度*********/
		if (sort[1]==2||sort[1]==4)
		{
			n_liquid=1;
			while (n_liquid<=n_liquid_total)
			{
				jump=0;
				/*液塞左界面净流入质量（正为流入，负为流出）=左侧气泡质量变化的一半*/
				dm_left[n_liquid]=-Tr_v[n_liquid]/2*dt;		
				dl_left[n_liquid]=dm_left[n_liquid]/A_l/ro_l;	//液塞左界面移动距离
				/*液塞右界面净流入质量（正为流入，负为流出）=右侧气泡质量变化的一半*/
				if (n_liquid==n_liquid_total)
					{dm_right[n_liquid]=-Tr_v[1]/2*dt;}
				else
					{dm_right[n_liquid]=-Tr_v[n_liquid+1]/2*dt;}
				dl_right[n_liquid]=dm_right[n_liquid]/A_l/ro_l;    //液塞右界面移动距离
				/*左界面控制体的重力系数，阻力系数，液体长度*/
				n=n_trans[2*n_liquid-1];	//定位到液塞左界面
				G_force[n_liquid]=0;	//重力初始化为0
				l[n_liquid]=0;
				C_ll=0;
				K=0;
				G_force[n_liquid]=G_force[n_liquid]+ro_l*g*c_gravity[n]*x_lright[n][1];	//将每个控制体的重力累加，x_lright[n][1]
				l[n_liquid]=l[n_liquid]+x_lright[n][1];	//将每个控制体的液体长度累加，x_lright[n][1]
				Re_lright[n]=d*ro_l*fabs(v_lright[n][1])/mu_right[n];	//控制体雷诺数，v_lright[n][1]
				if (Re_lright[n]==0)	//确定阻力系数，默认为正值
					{C_l[n]=0;}
				else if (Re_lright[n]<=1180)
					{C_l[n]=16/Re_lright[n];}
				else
					{C_l[n]=0.078*pow(Re_lright[n],-0.2);}
				if((Re_lright[n]>0)&&(n%210>=101)&&(n%210<=210))	//弯头处的阻力损失
					{K=1000/Re_lright[n]+0.1*(1+4/pow(2.3/25.4,0.3));
					C_l[n]=C_l[n]+K/4*2.3/110;}
				if (v_lright[n_trans[2*n_liquid-1]][1]>0)	//左界面速度大于0，即向正方向运动
					{C_l[n]=-C_l[n];}	//阻力系数为负，即阻力指向负方向
				C_ll=C_ll+C_l[n]*x_lright[n][1];	//控制体阻力系数乘以液体长度累加到总阻力系数
				/*左界面以右的所有控制体，重力系数，阻力系数，液塞长度*/
				n=n+1;
				while (n<=n_trans[2*n_liquid])
				{
					G_force[n_liquid]=G_force[n_liquid]+ro_l*g*c_gravity[n]*x_lleft[n][1];	//将每个控制体的重力累加，x_lleft[n][1]         
					l[n_liquid]=l[n_liquid]+x_lleft[n][1];	//将每个控制体的液体长度累加，x_lleft[n][1]
					Re_lleft[n]=d*ro_l*fabs(v_lleft[n][1])/mu_left[n];	//控制体雷诺数，v_lleft[n][1]
					if (Re_lleft[n]==0)
						{C_l[n]=0;}
					else if (Re_lleft[n]<=1180)
						{C_l[n]=16/Re_lleft[n];}
					else
						{C_l[n]=0.078*pow(Re_lleft[n],-0.2);}
					if((Re_lleft[n]>0)&&(n%210>=101)&&(n%210<=210))	//弯头处的阻力损失
						{K=1000/Re_lleft[n]+0.1*(1+4/pow(2.3/25.4,0.3));
						C_l[n]=C_l[n]+K/4*2.3/110;}
					if (v_lleft[n_trans[2*n_liquid-1]][1]>0)
						{C_l[n]=-C_l[n];}
					C_ll=C_ll+C_l[n]*x_lleft[n][1];
					n=n+1;
				}
				C_ll=C_ll/l[n_liquid];	//经过控制体液体长度加权平均的总阻力系数

				/**液塞动量方程，两边都除以了液塞截面积A_l，得到左界面控制体的速度**/
				v_lright[n_trans[2*n_liquid-1]][2]=1/(l[n_liquid]+dl_left[n_liquid]+dl_right[n_liquid]-dl_vl_left[n_liquid]-dl_vl_right[n_liquid])
						*(v_lright[n_trans[2*n_liquid-1]][1]*l[n_liquid]
						+dt/ro_l*(G_force[n_liquid]+P_v[n_trans[2*n_liquid-1]][1]-P_v[n_trans[2*n_liquid]][1]+2*C_ll*ro_l/d*l[n_liquid]*pow(v_lright[n_trans[2*n_liquid-1]][1],2)));

				/*将计算出左界面的v_l值赋给所有在同一液塞的控制体*/
				n=n_trans[2*n_liquid-1]+1;	//从左界面右边的控制体开始
				while (n<=n_trans[2*n_liquid]-1)
				{
					v_lright[n][2]=v_lright[n_trans[2*n_liquid-1]][2];
					v_lleft[n][2]=v_lright[n_trans[2*n_liquid-1]][2]; 
					n=n+1;
				}
				n=n_trans[2*n_liquid];	//右界面单独考虑，它不需要v_lright[n][2]
				v_lleft[n][2]=v_lright[n_trans[2*n_liquid-1]][2];
		
				/*下面判断是否有位置的偏移（液塞的起点或终点控制体位置改变）*/
				/*液塞的左界面*/
				n=n_trans[2*n_liquid-1];
				x_lright[n][2]=x_lright[n][1]-v_lright[n][1]*dt+dl_left[n_liquid]-dl_vl_left[n_liquid];	//左界面新的液体长度
				if (sort[n]!=8)	//不是液-气-液
					{x_lleft[n][2]=x_lleft[n][1];}	//左边液体长度保持不变
				/*若超出范围（向右或左跨越一个以上的控制体），跳出循环while (n_liquid<=n_liquid_total)*/
				if (x_lright[n][2]<=-dx||x_lright[n][2]>=2*dx)	
				{
					codenumber=1;	//记录异常
					break;	//跳出循环
				}
				/*正常情况，无偏移，左界面仍在原控制体内*/
				if (x_lright[n][2]>0&&x_lright[n][2]<dx)	
				{
					/*左界面液体的能量方程*/
					T_lright[n][2]=1/x_lright[n][2]*(x_lright[n][1]*T_lright[n][1]+dt/(c_pl_right[n]*ro_l*A_l)
						*(-v_lright[n][1]*A_l*ro_l*c_pl_right[n]*(x_lleft[n+1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[n+1][1])/(x_lright[n][1]+x_lleft[n+1][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]
						-2*A_l*(T_lright[n][1]-T_lleft[n+1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[n+1][1]/lamt_left[n+1])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));	//液体储能=流动带来的+壁面传热+轴向导热+蒸发冷凝-沸腾
					if (T_lright[n][2]>=T_max)
					{
						T_lright[n][2]=T_max;
					}
					if (T_lright[n][2]<=T_min)
					{
						T_lright[n][2]=T_min;
					}
					x_v[n][2]=dx-x_lright[n][2]-x_lleft[n][2];	//左界面新的气体长度
				}
				/*左界面向右移动一个控制体，将n、n+1合并为一个控制体*/
				else if (x_lright[n][2]<=0)
				{
					x_lright[n+1][2]=dx+x_lright[n][2];	//左界面新的液体长度
                    x_v[n+1][2]=dx-x_lright[n+1][2];	//左界面新的气体长度
                    x_lleft[n+1][2]=0;
                    T_lleft[n+1][2]=0;
                    v_lleft[n+1][2]=0;
                    sort[n+1]=4;	//4表示气-液
                    jump=1;	//记录左界面右移1个控制体
                      
                    /*左界面液体的能量方程*/
					T_lright[n+1][2]=1/(x_lright[n+1][2]*c_pl_right[n+1])*(x_lright[n][1]*T_lright[n][1]*c_pl_right[n]+x_lleft[n+1][1]*T_lleft[n+1][1]*c_pl_left[n+1]
						+dt/(A_l*ro_l)*(-v_lleft[n+1][1]*A_l*ro_l*c_pl_right[n+1]*(x_lleft[n+2][1]*T_lright[n+1][1]+x_lright[n+1][1]*T_lleft[n+2][1])/(x_lright[n+1][1]+x_lleft[n+2][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]+h_lleft[n+1]*(T_w[n+1][1]-T_lleft[n+1][1])*pi*d*x_lleft[n+1][1]
						-2*A_l*(T_lright[n+1][1]-T_lleft[n+2][1])/(x_lright[n+1][1]/lamt_right[n+1]+x_lleft[n+2][1]/lamt_left[n+2])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));	//液体储能=流动带来的+壁面传热+轴向导热+蒸发冷凝-沸腾，其中流动项和轴向导热项未考虑x_lright[n][1]???
					if (T_lright[n+1][2]>=T_max)
					{
						T_lright[n+1][2]=T_max;
					}
					if (T_lright[n+1][2]<=T_min)
					{
						T_lright[n+1][2]=T_min;
					}

					if (sort[n]!=8)	//不是液-气-液，即气-液
					{
						x_lright[n][2]=0;
                        x_lleft[n][2]=0;
                        x_v[n][2]=dx;
                        T_lleft[n][2]=0;
                        T_lright[n][2]=0;
                        v_lleft[n][2]=0;
                        v_lright[n][2]=0;
                        sort[n]=2;	//由于左界面移动变成了气
					}
					else  //n为液-气-液
					{
						x_lright[n][2]=0;
                        x_v[n][2]=dx-x_lleft[n][2];
                        sort[n]=3;	//由于左界面移动变成了液-气
                        T_lright[n][2]=0;
                        v_lright[n][2]=0;
					}
				}
				/*x_lright[n][2]>=dx，左界面向左移动一个控制体，将n、n-1合并为一个控制体*/
				else 
				{
					/*左界面液体的能量方程*/
					T_lright[n][2]=1/x_lright[n][2]*(x_lright[n][1]*T_lright[n][1]+dt/(c_pl_right[n]*ro_l*A_l)
						*(-v_lright[n][1]*A_l*ro_l*c_pl_right[n]*(x_lleft[n+1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[n+1][1])/(x_lright[n][1]+x_lleft[n+1][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]
						-2*A_l*(T_lright[n][1]-T_lleft[n+1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[n+1][1]/lamt_left[n+1])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));	//与左界面在原控制体内的方程一样
					if (T_lright[n][2]>=T_max)
					{
						T_lright[n][2]=T_max;
					}
					if (T_lright[n][2]<=T_min)
					{
						T_lright[n][2]=T_min;
					}

					if (n==1)	//左界面原来在第1个控制体，现在移动到最后一个控制体
					{
						if (sort[length]==3)	//最后一个控制体为液-气
						{
							T_lright[length][2]=T_lright[n][2];
                            v_lright[length][2]=v_lright[n][2];
                            x_lright[length][2]=x_lright[n][2]-dx;
                            x_v[length][2]=dx-x_lright[length][2]-x_lleft[length][2];
                            sort[length]=8;	//由于左界面移动变成液-气-液
						}
						else  //最后一个控制体为气
						{
							x_lright[length][2]=x_lright[n][2]-dx;
                            x_lleft[length][2]=0;
                            //x_lright[length][2]=x_lright[n][2]-dx;
                            T_lright[length][2]=T_lright[n][2];
                            T_lleft[length][2]=0;
                            v_lright[length][2]=v_lright[n][2];
                            v_lleft[length][2]=0;
                            sort[length]=4;	//由于左界面移动变成气-液
						}
						cross_negative=1;	//表示液塞左界面负向跨过原点，气泡序号发生变化
					}
					else  //左界面原来不在第1个控制体
					{
						if (sort[n-1]==3)	//液-气
						{
							T_lright[n-1][2]=T_lright[n][2];
                            v_lright[n-1][2]=v_lright[n][2];
                            x_lright[n-1][2]=x_lright[n][2]-dx;
                            x_v[n-1][2]=dx-x_lright[n-1][2]-x_lleft[n-1][2];
                            sort[n-1]=8;	//变成液-气-液
						}
						else  //气
						{
							x_lright[n-1][2]=x_lright[n][2]-dx;
                            x_lleft[n-1][2]=0;
                            //x_lright[n-1][2]=x_lright[n][2]-dx;
                            T_lright[n-1][2]=T_lright[n][2];
                            T_lleft[n-1][2]=0;
                            v_lright[n-1][2]=v_lright[n][2];
                            v_lleft[n-1][2]=0;
                            sort[n-1]=4;	//变成气-液
						}
					}

					x_lright[n][2]=dx;
                    x_lleft[n][2]=dx;
                    v_lleft[n][2]=v_lright[n][2];
                    T_lleft[n][2]=T_lright[n][2];
                    x_v[n][2]=0;
                    ro_v[n][2]=0;
                    T_v[n][2]=0;
					T_sat[n][2] = 0;
                    P_v[n][2]=0;
                    sort[n]=1;	//原左界面控制体变成液

				}

				/*液塞的中间部分*/
				if (jump==1)	//左界面向右移动一个控制体
				{
					n=n_trans[2*n_liquid-1]+2;	//定位到新左界面右侧的控制体
					while (n<=n_trans[2*n_liquid]-1)
					{
						/*液控制体的能量方程*/
						T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
							*((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
							-(x_lleft[n+1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[n+1][1])/(x_lright[n][1]+x_lleft[n+1][1]))
							+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
							-2*A_l*(T_lright[n][1]-T_lleft[n+1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[n+1][1]/lamt_left[n+1]));	//液体储存的能量=流动带来的+壁面传热+轴向导热
						if (T_lleft[n][2]>=T_max)
						{
							T_lleft[n][2]=T_max;
						}
						if (T_lleft[n][2]<=T_min)
						{
							T_lleft[n][2]=T_min;
						}
						T_lright[n][2]=T_lleft[n][2];
						x_lright[n][2]=dx;
						x_lleft[n][2]=dx;
						x_v[n][2]=0;
						n=n+1;
					}
				}
				else  //左界面在原控制体或向左移动一个控制体
				{
					n=n_trans[2*n_liquid-1]+1;	//定位到左界面右侧的控制体
					while (n<=n_trans[2*n_liquid]-1)
					{
						/*液控制体的能量方程*/
						T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
							*((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
							-(x_lleft[n+1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[n+1][1])/(x_lright[n][1]+x_lleft[n+1][1]))
							+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
							-2*A_l*(T_lright[n][1]-T_lleft[n+1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[n+1][1]/lamt_left[n+1]));	//同上
						if (T_lleft[n][2]>=T_max)
						{
							T_lleft[n][2]=T_max;
						}
						if (T_lleft[n][2]<=T_min)
						{
							T_lleft[n][2]=T_min;
						}
						T_lright[n][2]=T_lleft[n][2];
						x_lright[n][2]=dx;
						x_lleft[n][2]=dx;
						x_v[n][2]=0;
						n=n+1;
					}
				}
				jump=0;	//重置？？？为什么放这里

				/*液塞的右界面*/
				n=n_trans[2*n_liquid];
				x_lleft[n][2]=x_lleft[n][1]+v_lleft[n][1]*dt+dl_right[n_liquid]-dl_vl_right[n_liquid];	//右界面控制体液体长度
				if (sort[n]!=8)	//不是液-气-液，即只能是液-气
					{x_lright[n][2]=x_lright[n][1];}
				/*若超出范围（向右或左跨越一个以上的控制体），跳出循环while (n_liquid<=n_liquid_total)*/
				if (x_lleft[n][2]<=-dx||x_lleft[n][2]>=2*dx)
				{
					codenumber=1;	//记录异常
					break;	//跳出循环
				}
				/*正常情况，无偏移，右界面仍在原控制体内*/
				if (x_lleft[n][2]>0&&x_lleft[n][2]<dx)
				{
					/*液塞右界面的能量方程*/
					T_lleft[n][2]=1/x_lleft[n][2]*(x_lleft[n][1]*T_lleft[n][1]+dt/(c_pl_left[n]*ro_l*A_l)
						*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]*(x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lleft[n][1]+x_lright[n-1][1])
						+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
						-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
						+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));	//液体储能=流动带来的+壁面传热+轴向导热+蒸发冷凝-沸腾
					if (T_lleft[n][2]>=T_max)
					{
						T_lleft[n][2]=T_max;
					}
					if (T_lleft[n][2]<=T_min)
					{
						T_lleft[n][2]=T_min;
					}
					x_v[n][2]=dx-x_lleft[n][2]-x_lright[n][2];
				}
				/*右界面向右移动一个控制体，n、n+1合并为同一控制体*/
				else if (x_lleft[n][2]>=dx&&x_lleft[n][2]<2*dx)		
				{
					if (n==length)	//特殊情况，右界面从最后一个控制体移动到第一个控制体
					{
						/*右界面的能量方程*/
                        T_lleft[n][2]=1/x_lleft[n][2]*(x_lleft[n][1]*T_lleft[n][1]+dt/(c_pl_left[n]*ro_l*A_l)
							*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]*(x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lleft[n][1]+x_lright[n-1][1])
							+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
							-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
							+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));	//同上，全部没有考虑n+1控制体？？？
						if (T_lleft[n][2]>=T_max)
						{
							T_lleft[n][2]=T_max;
						}
						if (T_lleft[n][2]<=T_min)
						{
							T_lleft[n][2]=T_min;
						}

						if (sort[1]!=4)	//第1个控制体不为气-液，则只能为气
						{
							sort[1]=3;	//变成液-气
                            x_lleft[1][2]=x_lleft[n][2]-dx;
                            x_lright[1][2]=0;
                            x_v[1][2]=dx-x_lleft[1][2]-x_lright[1][2];
                            v_lleft[1][2]=v_lleft[n][2];
                            v_lright[1][2]=0;
                            T_lleft[1][2]=T_lleft[n][2];
                            T_lright[1][2]=0;
						}
						else  //第1个控制体为气-液
						{
							sort[1]=8;	//变成液-气-液
                            x_lleft[1][2]=x_lleft[n][2]-dx;
                            x_v[1][2]=dx-x_lleft[1][2]-x_lright[1][2];
                            v_lleft[1][2]=v_lleft[n][2];
                            T_lleft[1][2]=T_lleft[n][2];
						}
					}
					else  //右界面不在最后一个控制体
					{
						/*右界面的能量方程*/
                        T_lleft[n][2]=1/x_lleft[n][2]*(x_lleft[n][1]*T_lleft[n][1]+dt/(c_pl_left[n]*ro_l*A_l)
							*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]*(x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lleft[n][1]+x_lright[n-1][1])
							+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
							-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
							+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));	//同上上，全部没有考虑n+1控制体？？？
						if (T_lleft[n][2]>=T_max)
						{
							T_lleft[n][2]=T_max;
						}
						if (T_lleft[n][2]<=T_min)
						{
							T_lleft[n][2]=T_min;
						}

						if (sort[n+1]!=4)	//不是气-液，只能是气
						{
							sort[n+1]=3;	//变成液-气
                            x_lleft[n+1][2]=x_lleft[n][2]-dx;
                            x_lright[n+1][2]=0;
                            x_v[n+1][2]=dx-x_lleft[n+1][2]-x_lright[n+1][2];
                            v_lleft[n+1][2]=v_lleft[n][2];
                            v_lright[n+1][2]=0;
                            T_lleft[n+1][2]=T_lleft[n][2];
                            T_lright[n+1][2]=0;
						}
						else  //气-液
						{
							sort[n+1]=8;	//变成液-气-液
                            x_lleft[n+1][2]=x_lleft[n][2]-dx;
                            x_v[n+1][2]=dx-x_lleft[n+1][2]-x_lright[n+1][2];
                            v_lleft[n+1][2]=v_lleft[n][2];
                            T_lleft[n+1][2]=T_lleft[n][2];
						}
					}

					sort[n]=1;	//原右界面控制体变成液
                    x_lleft[n][2]=dx;
                    x_lright[n][2]=x_lleft[n][2];
                    v_lright[n][2]=v_lleft[n][2];
                    T_lright[n][2]=T_lleft[n][2];
                    x_v[n][2]=0;
                    ro_v[n][2]=0;
                    T_v[n][2]=0;
					T_sat[n][2] = 0;
                    P_v[n][2]=0;
				}
				/*右界面向左移动一个控制体，n、n-1合并为同一控制体*/
				else		
				{
					x_lleft[n-1][2]=dx+x_lleft[n][2];
                    x_lright[n-1][2]=0;
                    x_v[n-1][2]=dx-x_lleft[n-1][2];
                    T_lright[n-1][2]=0;
                    v_lright[n-1][2]=0;
                    sort[n-1]=3;	//变成液-气

					/*右界面的能量方程*/
                    T_lleft[n-1][2]=1/(x_lleft[n-1][2]*c_pl_left[n-1])*(x_lleft[n][1]*T_lleft[n][1]*c_pl_left[n]+x_lright[n-1][1]*T_lright[n-1][1]*c_pl_right[n-1]
						+dt/(A_l*ro_l)*(v_lright[n-1][1]*A_l*ro_l*c_pl_right[n-1]*(x_lright[n-2][1]*T_lleft[n-1][1]+x_lleft[n-1][1]*T_lright[n-2][1])/(x_lleft[n-1][1]+x_lright[n-2][1])
						+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]+h_lright[n-1]*(T_w[n-1][1]-T_lright[n-1][1])*pi*d*x_lright[n-1][1]
						-2*A_l*(T_lleft[n-1][1]-T_lright[n-2][1])/(x_lleft[n-1][1]/lamt_left[n-1]+x_lright[n-2][1]/lamt_right[n-2])
						+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));	//考虑了n-1和n的能量
					if (T_lleft[n-1][2]>=T_max)
					{
						T_lleft[n-1][2]=T_max;
					}
					if (T_lleft[n-1][2]<=T_min)
					{
						T_lleft[n-1][2]=T_min;
					}

					if (sort[n]==8)	//原右界面控制体为液-气-液
					{
						sort[n]=4;	//变成气-液
                        T_lleft[n][2]=0;
                        v_lleft[n][2]=0;
                        x_lleft[n][2]=0;
                        x_v[n][2]=dx-x_lright[n][2];
					}
					else  //原右界面控制体为液-气
					{
						sort[n]=2;	//变成气
                        T_lleft[n][2]=0;
                        T_lright[n][2]=0;
                        v_lleft[n][2]=0;
                        v_lright[n][2]=0;
                        x_lleft[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx;
					}
				}
				/*液塞从左界面到右界面讨论完毕，考虑下一个液塞*/
				n_liquid=n_liquid+1;
			}

			/*若超出范围，跳出外层循环while (i<=nn)*/
			if (codenumber==1)
				{break;}
		}

		/******对于起始点状态的分类讨论：第二部分：液、液-气、液-气-液，此时气泡序号领先于液塞序号******/
		else if (sort[1]==1||sort[1]==3||sort[1]==8)
		{
			/*确定各液塞的质量变化和重力作用、液塞总长度*/

			/***********************先对n_liquid==1做讨论，第1个液塞跨过了坐标原点************************/
			n_liquid=1;
			jump=0;
			dm_left[n_liquid]=-Tr_v[n_liquid_total]/2*dt;	//最后一个气泡相变质量的一半
			dl_left[n_liquid]=dm_left[n_liquid]/A_l/ro_l;	//左边界移动距离
			dm_right[n_liquid]=-Tr_v[n_liquid]/2*dt;	//第1个气泡相变质量的一半
			dl_right[n_liquid]=dm_right[n_liquid]/A_l/ro_l;	//右边界移动距离
        
			n=n_trans[2*n_liquid_total];	//定位到最后一个气液界面，即第1个液塞的左界面
			G_force[n_liquid]=0;
			l[n_liquid]=0;
			C_ll=0;
			K=0;
        
			G_force[n_liquid]=G_force[n_liquid]+ro_l*g*c_gravity[n]*x_lright[n][1];	//各控制体的重力累加
			l[n_liquid]=l[n_liquid]+x_lright[n][1];
			Re_lright[n]=d*ro_l*fabs(v_lright[n][1])/mu_right[n];
			if (Re_lright[n]==0)
				{C_l[n]=0;}
			else if (Re_lright[n]<=1180)
				{C_l[n]=16/Re_lright[n];}
			else 
				{C_l[n]=0.078*pow(Re_lright[n],-0.2);}
			if(Re_lright[n]>0&&(n%210>=101)&&(n%210<=210))	//弯头处的阻力损失
				{K=1000/Re_lright[n]+0.1*(1+4/pow(2.3/25.4,0.3));
				C_l[n]=C_l[n]+K/4*2.3/110;}
			if (v_lright[n_trans[2*n_liquid-1]][1]>0)	//速度向右则阻力向左
				{C_l[n]=-C_l[n];}
			C_ll=C_ll+C_l[n]*x_lright[n][1];	//液塞的阻力系数为各控制体阻力系数的加权平均

			n=n+1;
			/*从左界面到坐标原点这部分液塞*/
			while (n<=length)
			{
				G_force[n_liquid]=G_force[n_liquid]+ro_l*g*c_gravity[n]*x_lleft[n][1];	//重力作用
				l[n_liquid]=l[n_liquid]+x_lleft[n][1];	//液塞长度
				Re_lleft[n]=d*ro_l*fabs(v_lleft[n][1])/mu_left[n];
				if (Re_lleft[n]==0)
					{C_l[n]=0;}
				else if (Re_lleft[n]<=1180)
					{C_l[n]=16/Re_lleft[n];}
				else 
					{C_l[n]=0.078*pow(Re_lleft[n],-0.2);}
				if(Re_lleft[n]>0&&(n%210>=101)&&(n%210<=210))	//弯头处的阻力损失
					{K=1000/Re_lleft[n]+0.1*(1+4/pow(2.3/25.4,0.3));
					C_l[n]=C_l[n]+K/4*2.3/110;}
				if (v_lleft[n_trans[2*n_liquid-1]][1]>0)
					{C_l[n]=-C_l[n];}
				C_ll=C_ll+C_l[n]*x_lleft[n][1];
				n=n+1;
			}
			n=1;
			/*从坐标原点到右界面这部分液塞*/
			while (n<=n_trans[2*n_liquid-1])
			{
				G_force[n_liquid]=G_force[n_liquid]+ro_l*g*c_gravity[n]*x_lleft[n][1];
				l[n_liquid]=l[n_liquid]+x_lleft[n][1];                      
				Re_lleft[n]=d*ro_l*fabs(v_lleft[n][1])/mu_left[n];
				if (Re_lleft[n]==0)
					{C_l[n]=0;}
				else if (Re_lleft[n]<=1180)
					{C_l[n]=16/Re_lleft[n];}
				else 
					{C_l[n]=0.078*pow(Re_lleft[n],-0.2);}
				if(Re_lleft[n]>0&&(n%210>=101)&&(n%210<=210))	//弯头处的阻力损失
					{K=1000/Re_lleft[n]+0.1*(1+4/pow(2.3/25.4,0.3));
					C_l[n]=C_l[n]+K/4*2.3/110;}
				if (v_lleft[n_trans[2*n_liquid-1]][1]>0)
					{C_l[n]=-C_l[n];}
				C_ll=C_ll+C_l[n]*x_lleft[n][1];
				n=n+1;
			}

			C_ll=C_ll/l[n_liquid];	//之前加权累加，此处平均

			/**动量方程，两边都除以了液塞截面积A_l，得到右界面控制体的速度**/
			v_lleft[n_trans[2*n_liquid-1]][2]=1/(l[n_liquid]+dl_left[n_liquid]+dl_right[n_liquid]-dl_vl_left[n_liquid]-dl_vl_right[n_liquid])
				*(v_lleft[n_trans[2*n_liquid-1]][1]*l[n_liquid]+dt/ro_l*(G_force[n_liquid]+P_v[n_trans[2*n_liquid_total]][1]-P_v[n_trans[2*n_liquid-1]][1]
				+2*C_ll*ro_l/d*l[n_liquid]*pow(v_lleft[n_trans[2*n_liquid-1]][1],2)));

			/*将计算出的v_l值赋给所有同一液塞的控制体*/
			n=n_trans[2*n_liquid_total];	//定位到左界面
			while (n<=length)	//从左界面到坐标原点
			{
				if (n!=n_trans[2*n_liquid_total])
				{
					v_lleft[n][2]=v_lleft[n_trans[2*n_liquid-1]][2]; 
					v_lright[n][2]=v_lleft[n_trans[2*n_liquid-1]][2];
				}
				else  //左界面没有左边部分液塞
				{
					v_lright[n][2]=v_lleft[n_trans[2*n_liquid-1]][2];
				}
				n=n+1;
			}

			n=1;	//定位到第1个控制体
			while (n<=n_trans[2*n_liquid-1])	//从坐标原点到右界面
			{
				if (n!=n_trans[2*n_liquid-1])
				{
					v_lleft[n][2]=v_lleft[n_trans[2*n_liquid-1]][2];
					v_lright[n][2]=v_lleft[n_trans[2*n_liquid-1]][2];
				}
				else  //右界面没有右边部分液塞
				{
					v_lleft[n][2]=v_lleft[n_trans[2*n_liquid-1]][2];
				}
				n=n+1;
			}

			/*下面判断是否有位置的偏移（液塞的起点或终点控制体位置改变）*/

			/*左界面*/
			n=n_trans[2*n_liquid_total];
			x_lright[n][2]=x_lright[n][1]-v_lright[n][1]*dt+dl_left[n_liquid]-dl_vl_left[n_liquid];	//左界面控制体新的液体长度
			if (sort[n]!=8)	//如果是液-气-液，则左边部分液体长度保持不变
				{x_lleft[n][2]=x_lleft[n][1];}
			/*若左界面移动超过一个控制体，跳出循环while (i<=nn)*/
			if (x_lright[n][2]<=-dx||x_lright[n][2]>=2*dx)	
			{
				codenumber=2;
				break;
			}

			/*左界面正常情况，无偏移*/
			if (x_lright[n][2]>0&&x_lright[n][2]<dx)			
			{
				/*左界面液体的能量方程*/
				if (n!=length)
				{
					T_lright[n][2]=1/x_lright[n][2]*(x_lright[n][1]*T_lright[n][1]+dt/(c_pl_right[n]*ro_l*A_l)
						*(-v_lright[n][1]*A_l*ro_l*c_pl_right[n]*(x_lleft[n+1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[n+1][1])/(x_lright[n][1]+x_lleft[n+1][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]
						-2*A_l*(T_lright[n][1]-T_lleft[n+1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[n+1][1]/lamt_left[n+1])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));
					if (T_lright[n][2]>=T_max)
					{
						T_lright[n][2]=T_max;
					}
					if (T_lright[n][2]<=T_min)
					{
						T_lright[n][2]=T_min;
					}
					x_v[n][2]=dx-x_lright[n][2]-x_lleft[n][2];
				}
				else  //如果左界面在最后一个控制体，则n+1=1
				{
					T_lright[n][2]=1/x_lright[n][2]*(x_lright[n][1]*T_lright[n][1]+dt/(c_pl_right[n]*ro_l*A_l)
						*(-v_lright[n][1]*A_l*ro_l*c_pl_right[n]*(x_lleft[1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[1][1])/(x_lright[n][1]+x_lleft[1][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]
						-2*A_l*(T_lright[n][1]-T_lleft[1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[1][1]/lamt_left[1])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));
					if (T_lright[n][2]>=T_max)
					{
						T_lright[n][2]=T_max;
					}
					if (T_lright[n][2]<=T_min)
					{
						T_lright[n][2]=T_min;
					}
					x_v[n][2]=dx-x_lright[n][2]-x_lleft[n][2];
				}
			}
			/*左界面向右移动一个控制体，将n、n+1合并为一个控制体*/
			else if (x_lright[n][2]<=0)							
			{
				/*左界面原来在最后一个控制体，移动至第1个控制体*/
				if (n==length)
				{
					x_lright[1][2]=dx+x_lright[n][2];
                    x_lleft[1][2]=0;
                    x_v[1][2]=dx-x_lright[1][2]-x_lleft[1][2];
                    T_lleft[1][2]=0;
                    v_lleft[1][2]=0;
                    sort[1]=4;	//第1个控制体变成气-液
                    jump=1;
                    cross_positive=1;	//液塞左界面正向跨过原点，气泡序号加1

					/*左界面液体的能量方程*/
                    T_lright[1][2]=1/(x_lright[1][2]*c_pl_right[1])*(x_lright[n][1]*T_lright[n][1]*c_pl_right[n]+x_lleft[1][1]*T_lleft[1][1]*c_pl_left[1]
						+dt/(A_l*ro_l)*(-v_lright[1][1]*A_l*ro_l*c_pl_right[1]*(x_lleft[2][1]*T_lright[1][1]+x_lright[1][1]*T_lleft[2][1])/(x_lright[1][1]+x_lleft[2][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]+h_lright[1]*(T_w[1][1]-T_lright[1][1])*pi*d*x_lright[1][1]
						-2*A_l*(T_lright[1][1]-T_lleft[2][1])/(x_lright[1][1]/lamt_right[1]+x_lleft[2][1]/lamt_left[2])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));
                    if (T_lright[1][2]>=T_max)
					{
						T_lright[1][2]=T_max;
					}
					if (T_lright[1][2]<=T_min)
					{
						T_lright[1][2]=T_min;
					}

					if (sort[n]!=8)	//左界面控制体原本为气-液
					{
						T_lleft[n][2]=0;
                        T_lright[n][2]=0;
                        v_lleft[n][2]=0;
                        v_lright[n][2]=0;
                        x_lleft[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx;
                        sort[n]=2;	//界面移动后变成气
					}
					else  //左界面控制体原本为液-气-液
					{
						sort[n]=3;	//界面移动后变成液-气
                        T_lright[n][2]=0;
                        v_lright[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx-x_lleft[n][2];
					}
				}
				/*左界面原来在倒数第2个控制体，移动至最后1个控制体*/
				else if (n==(length-1))
				{
					x_lright[n+1][2]=dx+x_lright[n][2];
                    x_lleft[n+1][2]=0;
                    x_v[n+1][2]=dx-x_lright[n+1][2]-x_lleft[n+1][2];
                    T_lleft[n+1][2]=0;
                    v_lleft[n+1][2]=0;
                    sort[n+1]=4;	//最后一个控制体变成气-液
                    jump=1;

					/*左界面液体的能量方程*/
                    T_lright[n+1][2]=1/(x_lright[n+1][2]*c_pl_right[n+1])*(x_lright[n][1]*T_lright[n][1]*c_pl_right[n]+x_lright[n+1][1]*T_lright[n+1][1]*c_pl_right[n+1]
						+dt/(A_l*ro_l)*(-v_lright[n+1][1]*A_l*c_pl_right[n+1]*ro_l*(x_lleft[1][1]*T_lright[n+1][1]+x_lright[n+1][1]*T_lleft[1][1])/(x_lright[n+1][1]+x_lleft[1][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]+h_lright[n+1]*(T_w[n+1][1]-T_lright[n+1][1])*pi*d*x_lright[n+1][1]
						-2*A_l*(T_lright[n+1][1]-T_lleft[1][1])/(x_lright[n+1][1]/lamt_right[n+1]+x_lleft[1][1]/lamt_left[1])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));
                    if (T_lright[n+1][2]>=T_max)
					{
						T_lright[n+1][2]=T_max;
					}
					if (T_lright[n+1][2]<=T_min)
					{
						T_lright[n+1][2]=T_min;
					}

					if (sort[n]!=8)	//左界面控制体原本为气-液
					{
						T_lleft[n][2]=0;
                        T_lright[n][2]=0;
                        v_lleft[n][2]=0;
                        v_lright[n][2]=0;
                        x_lleft[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx;
                        sort[n]=2;	//现在变成气
					}
					else//左界面控制体原本为液-气-液
					{
						sort[n]=3;	//现在变成液-气
                        T_lright[n][2]=0;
                        v_lright[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx-x_lleft[n][2];
					}
				}
				/*左界面原来在倒数第3个控制体或者更左*/
				else
				{
					x_lright[n+1][2]=dx+x_lright[n][2];
                    x_lleft[n+1][2]=0;
                    x_v[n+1][2]=dx-x_lright[n+1][2]-x_lleft[n+1][2];
                    T_lleft[n+1][2]=0;
                    v_lleft[n+1][2]=0;
                    sort[n+1]=4;	//n+1变成现在的左界面控制体，气-液
                    jump=1;

					/*左界面液体的能量方程*/
                    T_lright[n+1][2]=1/(x_lright[n+1][2]*c_pl_right[n+1])*(x_lright[n][1]*T_lright[n][1]*c_pl_right[n]+x_lright[n+1][1]*T_lright[n+1][1]*c_pl_right[n+1]
						+dt/(A_l*ro_l)*(-v_lright[n+1][1]*A_l*ro_l*c_pl_right[n+1]*(x_lleft[n+2][1]*T_lright[n+1][1]+x_lright[n+1][1]*T_lleft[n+2][1])/(x_lright[n+1][1]+x_lleft[n+2][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]+h_lright[n+1]*(T_w[n+1][1]-T_lright[n+1][1])*pi*d*x_lright[n+1][1]
						-2*A_l*(T_lright[n+1][1]-T_lleft[n+2][1])/(x_lright[n+1][1]/lamt_right[n+1]+x_lleft[n+2][1]/lamt_left[n+2])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));
                    if (T_lright[n+1][2]>=T_max)
					{
						T_lright[n+1][2]=T_max;
					}
					if (T_lright[n+1][2]<=T_min)
					{
						T_lright[n+1][2]=T_min;
					}

					if (sort[n]!=8)	//原来为气-液
					{
						T_lleft[n][2]=0;
                        T_lright[n][2]=0;
                        v_lleft[n][2]=0;
                        v_lright[n][2]=0;
                        x_lleft[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx;
                        sort[n]=2;	//现变成气
					}
					else  //原来为液-气-液
					{
						sort[n]=3;	//现在为液-气
                        T_lright[n][2]=0;
                        v_lright[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx-x_lleft[n][2];
					}
				}
			}
			/*左界面向左移动一个控制体，将n、n-1合并为一个控制体*/
			else
			{
				/*左界面液体的能量方程*/
				/*左界面原来在最后一个控制体，移动至倒数第2个控制体*/
				if (n==length)
				{
					T_lright[n][2]=1/x_lright[n][2]*(x_lright[n][1]*T_lright[n][1]+dt/(c_pl_right[n]*ro_l*A_l)
						*(-v_lright[n][1]*A_l*ro_l*c_pl_right[n]*(x_lleft[1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[1][1])/(x_lright[n][1]+x_lleft[1][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]
						-2*A_l*(T_lright[n][1]-T_lleft[1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[1][1]/lamt_left[1])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));
				}
				/*左界面原来在倒数第2个控制体或者更左*/
				else
				{
					T_lright[n][2]=1/x_lright[n][2]*(x_lright[n][1]*T_lright[n][1]+dt/(c_pl_right[n]*ro_l*A_l)
						*(-v_lright[n][1]*A_l*ro_l*c_pl_right[n]*(x_lleft[n+1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[n+1][1])/(x_lright[n][1]+x_lleft[n+1][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]
						-2*A_l*(T_lright[n][1]-T_lleft[n+1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[n+1][1]/lamt_left[n+1])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));
				}
				if (T_lright[n][2]>=T_max)
				{
					T_lright[n][2]=T_max;
				}
				if (T_lright[n][2]<=T_min)
				{
					T_lright[n][2]=T_min;
				}

				if (sort[n-1]!=3)	//原左界面左边控制体为气
				{
					x_lright[n-1][2]=x_lright[n][2]-dx; 
                    x_lleft[n-1][2]=0;
                    x_v[n-1][2]=dx-x_lright[n-1][2]-x_lleft[n-1][2];
                    v_lright[n-1][2]=v_lright[n][2];
                    v_lleft[n-1][2]=0;
                    T_lright[n-1][2]=T_lright[n][2];
                    T_lleft[n-1][2]=0;
                    sort[n-1]=4;	//现变为气-液
                          
                    x_lright[n][2]=dx;
                    x_lleft[n][2]=dx;
                    x_v[n][2]=0;
                    v_lleft[n][2]=v_lright[n][2];
                    T_lleft[n][2]=T_lright[n][2];
                    ro_v[n][2]=0;
                    T_v[n][2]=0;
					T_sat[n][2] = 0;
                    P_v[n][2]=0;
                    sort[n]=1;	//原左界面变为液
				}
				else //原左界面左边控制体为液-气
				{
					x_lright[n-1][2]=x_lright[n][2]-dx; 
                    x_v[n-1][2]=dx-x_lright[n-1][2]-x_lleft[n-1][2];
                    v_lright[n-1][2]=v_lright[n][2];
                    T_lright[n-1][2]=T_lright[n][2];
                    sort[n-1]=8;	//现变为液-气-液

					x_lright[n][2]=dx;
                    x_lleft[n][2]=dx;
                    x_v[n][2]=0;
                    v_lleft[n][2]=v_lright[n][2];
                    T_lleft[n][2]=T_lright[n][2];
                    ro_v[n][2]=0;
                    T_v[n][2]=0;
					T_sat[n][2] = 0;
                    P_v[n][2]=0;
                    sort[n]=1;	//原左界面变为液
				}

			}

			/*液塞的中间部分*/
			if (jump==1)	//左界面向右移出控制体
			{
				if (n_trans[2*n_liquid_total]==length)	//左界面原来在最后一个控制体，现移至第1个控制体
				{
					n=2;
					while (n<=n_trans[2*n_liquid-1]-1)	//从第2个控制体到右界面左边控制体
					{
						/*液控制体的能量方程，两边同除以dx*/
						T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
							*((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
							-(x_lleft[n+1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[n+1][1])/(x_lleft[n][1]+x_lleft[n+1][1]))
							+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
							-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
							-2*A_l*(T_lleft[n][1]-T_lleft[n+1][1])/(x_lleft[n][1]/lamt_left[n]+x_lleft[n+1][1]/lamt_left[n+1]));
						if (T_lleft[n][2]>=T_max)
						{
							T_lleft[n][2]=T_max;
						}
						if (T_lleft[n][2]<=T_min)
						{
							T_lleft[n][2]=T_min;
						}
						T_lright[n][2]=T_lleft[n][2];
						x_lleft[n][2]=dx;
						x_lright[n][2]=dx;
						x_v[n][2]=0;
						n=n+1;
					}
				}
				else if (n_trans[2*n_liquid_total]==(length-1))	//左界面原来在倒数第2个控制体，现移动至最后一个控制体
				{
					n=1;
					while (n<=n_trans[2*n_liquid-1]-1)	//从第1个控制体到右界面左边控制体
					{
						if (n==1)	//第1个控制体单独处理
						{
							/*液控制体的能量方程*/
							T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
								*((x_lright[length][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[length][1])/(x_lright[length][1]+x_lleft[n][1])
								-(x_lleft[n+1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[n+1][1])/(x_lleft[n][1]+x_lleft[n+1][1]))
								+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
								-2*A_l*(T_lleft[n][1]-T_lright[length][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[length][1]/lamt_right[length])
								-2*A_l*(T_lleft[n][1]-T_lleft[n+1][1])/(x_lleft[n][1]/lamt_left[n]+x_lleft[n+1][1]/lamt_left[n+1]));
							if (T_lleft[n][2]>=T_max)
							{
								T_lleft[n][2]=T_max;
							}
							if (T_lleft[n][2]<=T_min)
							{
								T_lleft[n][2]=T_min;
							}
							T_lright[n][2]=T_lleft[n][2];
							x_lleft[n][2]=dx;
							x_lright[n][2]=dx;
							x_v[n][2]=0;
						}
						else  //从第2个控制体到右界面左边控制体
						{
							/*液控制体的能量方程*/
							T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]*
								((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
								-(x_lleft[n+1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[n+1][1])/(x_lleft[n][1]+x_lleft[n+1][1]))
								+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
								-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
								-2*A_l*(T_lleft[n][1]-T_lleft[n+1][1])/(x_lleft[n][1]/lamt_left[n]+x_lleft[n+1][1]/lamt_left[n+1]));
							if (T_lleft[n][2]>=T_max)
							{
								T_lleft[n][2]=T_max;
							}
							if (T_lleft[n][2]<=T_min)
							{
								T_lleft[n][2]=T_min;
							}
							T_lright[n][2]=T_lleft[n][2];
							x_lleft[n][2]=dx;
							x_lright[n][2]=dx;
							x_v[n][2]=0;
						}
						n=n+1;
					}
				}
				else  //左界面原来在倒数第3个控制体或者更左，将这个液塞以原点为界分为两段
				{
					n=n_trans[2*n_liquid_total]+2;	//新的左界面右边控制体，到最后一个控制体
					while (n<=length)
					{
						if (n==length)	//最后一个控制体单独处理
						{
							/*液塞的能量方程*/
							T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
								*((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
								-(x_lleft[1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[1][1])/(x_lleft[n][1]+x_lleft[1][1]))
								+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
								-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
								-2*A_l*(T_lleft[n][1]-T_lleft[1][1])/(x_lleft[n][1]/lamt_left[n]+x_lleft[1][1]/lamt_left[n+1]));
							if (T_lleft[n][2]>=T_max)
							{
								T_lleft[n][2]=T_max;
							}
							if (T_lleft[n][2]<=T_min)
							{
								T_lleft[n][2]=T_min;
							}
							T_lright[n][2]=T_lleft[n][2];
							x_lleft[n][2]=dx;
							x_lright[n][2]=dx;
							x_v[n][2]=0;
						}
						else  //新的左界面右边控制体，到倒数第2个控制体
						{
							/*液塞的能量方程*/
							T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
								*((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
								-(x_lleft[n+1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[n+1][1])/(x_lleft[n][1]+x_lleft[n+1][1]))
								+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
								-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
								-2*A_l*(T_lleft[n][1]-T_lleft[n+1][1])/(x_lleft[n][1]/lamt_left[n]+x_lleft[n+1][1]/lamt_left[n+1]));
							if (T_lleft[n][2]>=T_max)
							{
								T_lleft[n][2]=T_max;
							}
							if (T_lleft[n][2]<=T_min)
							{
								T_lleft[n][2]=T_min;
							}
							T_lright[n][2]=T_lleft[n][2];
							x_lleft[n][2]=dx;
							x_lright[n][2]=dx;
							x_v[n][2]=0;
						}
						n=n+1;
					}
					n=1;	//从第1个控制体到右界面
					if (n_trans[1]>=2)	//右界面在第2个控制体或者更右，否则第1个控制体就是右界面，在后面计算
					{
						while (n<=n_trans[1]-1)	//从第1个控制体到右界面左边控制体
						{
							if (n==1)	//第1个控制体单独处理
							{
								/*液塞的能量方程*/
								T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
									*((x_lright[length][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[length][1])/(x_lright[length][1]+x_lleft[n][1])
									-(x_lleft[n+1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[n+1][1])/(x_lleft[n][1]+x_lleft[n+1][1]))
									+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
									-2*A_l*(T_lleft[n][1]-T_lright[length][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[length][1]/lamt_right[length])
									-2*A_l*(T_lleft[n][1]-T_lleft[n+1][1])/(x_lleft[n][1]/lamt_left[n]+x_lleft[n+1][1]/lamt_left[n+1]));
								if (T_lleft[n][2]>=T_max)
								{
									T_lleft[n][2]=T_max;
								}
								if (T_lleft[n][2]<=T_min)
								{
									T_lleft[n][2]=T_min;
								}
								T_lright[n][2]=T_lleft[n][2];
								x_lleft[n][2]=dx;
								x_lright[n][2]=dx;
								x_v[n][2]=0;
							}
							else //从第2个控制体到右界面左边控制体
							{
								/*液塞的能量方程*/
								T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
									*((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
									-(x_lleft[n+1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[n+1][1])/(x_lleft[n][1]+x_lleft[n+1][1]))
									+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
									-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
									-2*A_l*(T_lleft[n][1]-T_lleft[n+1][1])/(x_lleft[n][1]/lamt_left[n]+x_lleft[n+1][1]/lamt_left[n+1]));
								if (T_lleft[n][2]>=T_max)
								{
									T_lleft[n][2]=T_max;
								}
								if (T_lleft[n][2]<=T_min)
								{
									T_lleft[n][2]=T_min;
								}
								T_lright[n][2]=T_lleft[n][2];
								x_lleft[n][2]=dx;
								x_lright[n][2]=dx;
								x_v[n][2]=0;
							}
							n=n+1;
						}
					}
				}
			}
			else  //左界面在原控制体或者向左移出控制体
			{
				if (n_trans[2*n_liquid_total]==length)	//左界面在最后一个控制体
				{
					n=1;
					while (n<=n_trans[2*n_liquid-1]-1)	//从第1个控制体到右界面左边控制体
					{
						if (n==1)	//第1个控制体单独处理
						{
							/*液塞的能量方程*/
							T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
								*((x_lright[length][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[length][1])/(x_lright[length][1]+x_lleft[n][1])
								-(x_lleft[n+1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[n+1][1])/(x_lleft[n][1]+x_lleft[n+1][1]))
								+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
								-2*A_l*(T_lleft[n][1]-T_lright[length][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[length][1]/lamt_right[length])
								-2*A_l*(T_lleft[n][1]-T_lleft[n+1][1])/(x_lleft[n][1]/lamt_left[n]+x_lleft[n+1][1]/lamt_left[n+1]));
							if (T_lleft[n][2]>=T_max)
							{
								T_lleft[n][2]=T_max;
							}
							if (T_lleft[n][2]<=T_min)
							{
								T_lleft[n][2]=T_min;
							}
							T_lright[n][2]=T_lleft[n][2];
							x_lleft[n][2]=dx;
							x_lright[n][2]=dx;
							x_v[n][2]=0;
						}
						else  //从第2个控制体到右界面左边控制体
						{
							/*液塞的能量方程*/
							T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
								*((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
								-(x_lleft[n+1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[n+1][1])/(x_lleft[n][1]+x_lleft[n+1][1]))
								+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
								-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
								-2*A_l*(T_lleft[n][1]-T_lleft[n+1][1])/(x_lleft[n][1]/lamt_left[n]+x_lleft[n+1][1]/lamt_left[n+1]));
							if (T_lleft[n][2]>=T_max)
							{
								T_lleft[n][2]=T_max;
							}
							if (T_lleft[n][2]<=T_min)
							{
								T_lleft[n][2]=T_min;
							}
							T_lright[n][2]=T_lleft[n][2];
							x_lleft[n][2]=dx;
							x_lright[n][2]=dx;
							x_v[n][2]=0;
						}
						n=n+1;
					}
				}
				else  //左界面不在最后一个控制体，将液塞以原点为界分段两段
				{
					n=n_trans[2*n_liquid_total]+1;
					while (n<=length)	//从左界面右边控制体到最后一个控制体
					{
						if (n==length)	//最后一个控制体单独处理
						{
							/*液塞的能量方程*/
							T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
								*((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
								-(x_lleft[1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[1][1])/(x_lleft[n][1]+x_lleft[1][1]))
								+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
								-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
								-2*A_l*(T_lleft[n][1]-T_lleft[1][1])/(x_lleft[n][1]/lamt_left[n]+x_lleft[1][1]/lamt_left[1]));
							if (T_lleft[n][2]>=T_max)
							{
								T_lleft[n][2]=T_max;
							}
							if (T_lleft[n][2]<=T_min)
							{
								T_lleft[n][2]=T_min;
							}
							T_lright[n][2]=T_lleft[n][2];
							x_lleft[n][2]=dx;
							x_lright[n][2]=dx;
							x_v[n][2]=0;
						}
						else  //从左界面右边控制体到倒数第2个控制体
						{
							/*液塞的能量方程*/
							T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
								*((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
								-(x_lleft[n+1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[n+1][1])/(x_lleft[n][1]+x_lleft[n+1][1]))
								+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
								-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
								-2*A_l*(T_lleft[n][1]-T_lleft[n+1][1])/(x_lleft[n][1]/lamt_left[n]+x_lleft[n+1][1]/lamt_left[n+1]));
							if (T_lleft[n][2]>=T_max)
							{
								T_lleft[n][2]=T_max;
							}
							if (T_lleft[n][2]<=T_min)
							{
								T_lleft[n][2]=T_min;
							}
							T_lright[n][2]=T_lleft[n][2];
							x_lleft[n][2]=dx;
							x_lright[n][2]=dx;
							x_v[n][2]=0;
						}
						n=n+1;
					}
					n=1;	//从第1个控制体到右界面左边控制体
					if (n_trans[1]>=2)	//右边界在第2个控制体或者更右边，否则第1个控制体就是右边界，在后面计算
					{
						while (n<=n_trans[1]-1)
						{
							if (n==1)	//第1个控制体单独处理
							{
								/*液塞的能量方程*/
								T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
									*((x_lright[length][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[length][1])/(x_lright[length][1]+x_lleft[n][1])
									-(x_lleft[n+1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[n+1][1])/(x_lleft[n][1]+x_lleft[n+1][1]))
									+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
									-2*A_l*(T_lleft[n][1]-T_lright[length][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[length][1]/lamt_right[length])
									-2*A_l*(T_lleft[n][1]/lamt_left[n]-T_lleft[n+1][1]/lamt_left[n+1])/(x_lleft[n][1]+x_lleft[n+1][1]));
								if (T_lleft[n][2]>=T_max)
								{
									T_lleft[n][2]=T_max;
								}
								if (T_lleft[n][2]<=T_min)
								{
									T_lleft[n][2]=T_min;
								}
								T_lright[n][2]=T_lleft[n][2];
								x_lleft[n][2]=dx;
								x_lright[n][2]=dx;
								x_v[n][2]=0;
							}
							else  //从第2个控制体到右界面左边控制体
							{
								/*液塞的能量方程*/
								T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
									*((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
									-(x_lleft[n+1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[n+1][1])/(x_lleft[n][1]+x_lleft[n+1][1]))
									+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
									-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
									-2*A_l*(T_lleft[n][1]/lamt_left[n]-T_lleft[n+1][1]/lamt_left[n+1])/(x_lleft[n][1]+x_lleft[n+1][1]));
								if (T_lleft[n][2]>=T_max)
								{
									T_lleft[n][2]=T_max;
								}
								if (T_lleft[n][2]<=T_min)
								{
									T_lleft[n][2]=T_min;
								}
								T_lright[n][2]=T_lleft[n][2];
								x_lleft[n][2]=dx;
								x_lright[n][2]=dx;
								x_v[n][2]=0;
							}
							n=n+1;
						}
					}
				}
			}

			jump=0;	//重置界面移动标记

			/*右界面*/
			n=n_trans[2*n_liquid-1];
			x_lleft[n][2]=x_lleft[n][1]+v_lleft[n][1]*dt+dl_right[n_liquid]-dl_vl_right[n_liquid];	//右界面控制体新的液体长度
			if (sort[n]!=8)	//右界面不是液-气-液，而是液-气
				{x_lright[n][2]=x_lright[n][1];}
			/*若右界面移动超过一个控制体，跳出循环while (i<=nn)*/
			if (x_lleft[n][2]<=-dx||x_lleft[n][2]>=2*dx)
			{
				codenumber=2;
				break;
			}
			/*右界面正常无偏移*/
			if (x_lleft[n][2]>0&&x_lleft[n][2]<dx)
			{
				if (n==1)//右界面在第1个控制体
				{
					/*右界面液塞的能量方程*/
					T_lleft[n][2]=1/x_lleft[n][2]*(x_lleft[n][1]*T_lleft[n][1]+dt/(c_pl_left[n]*ro_l*A_l)
						*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]*(x_lright[length][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[length][1])/(x_lleft[n][1]+x_lright[length][1])
						+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
						-2*A_l*(T_lleft[n][1]-T_lright[length][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[length][1]/lamt_right[length])
						+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));
					if (T_lleft[n][2]>=T_max)
					{
						T_lleft[n][2]=T_max;
					}
					if (T_lleft[n][2]<=T_min)
					{
						T_lleft[n][2]=T_min;
					}
					x_v[n][2]=dx-x_lleft[n][2]-x_lright[n][2];
				}
				else//右界面在第2个控制体或更右
				{
					/*右界面液塞的能量方程*/
					T_lleft[n][2]=1/x_lleft[n][2]*(x_lleft[n][1]*T_lleft[n][1]+dt/(c_pl_left[n]*ro_l*A_l)
						*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]*(x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lleft[n][1]+x_lright[n-1][1])
						+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
						-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
						+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));
					if (T_lleft[n][2]>=T_max)
					{
						T_lleft[n][2]=T_max;
					}
					if (T_lleft[n][2]<=T_min)
					{
						T_lleft[n][2]=T_min;
					}
					x_v[n][2]=dx-x_lleft[n][2]-x_lright[n][2];
				}
			}
			/*右界面向右移动一个控制体，n、n+1合并考虑*/
			else if (x_lleft[n][2]>=dx&&x_lleft[n][2]<2*dx)		
			{
				/*右界面液塞的能量方程*/
				if (n==1)//右界面原来在第1个控制体
				{
					T_lleft[n][2]=1/x_lleft[n][2]*(x_lleft[n][1]*T_lleft[n][1]+dt/(c_pl_left[n]*ro_l*A_l)
						*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]*(x_lright[length][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[length][1])/(x_lleft[n][1]+x_lright[length][1])
						+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
						-2*A_l*(T_lleft[n][1]-T_lright[length][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[length][1]/lamt_right[length])
						+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));
					if (T_lleft[n][2]>=T_max)
					{
						T_lleft[n][2]=T_max;
					}
					if (T_lleft[n][2]<=T_min)
					{
						T_lleft[n][2]=T_min;
					}
				}
				else//右界面原来在第2个控制体或更右
				{
					T_lleft[n][2]=1/x_lleft[n][2]*(x_lleft[n][1]*T_lleft[n][1]+dt/(c_pl_left[n]*ro_l*A_l)
						*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]*(x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lleft[n][1]+x_lright[n-1][1])
						+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
						-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
						+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));
					if (T_lleft[n][2]>=T_max)
					{
						T_lleft[n][2]=T_max;
					}
					if (T_lleft[n][2]<=T_min)
					{
						T_lleft[n][2]=T_min;
					}
				}

				if (sort[n+1]!=4)//原来为气
				{
					sort[n+1]=3;	//现在变成液-气
                    x_lleft[n+1][2]=x_lleft[n][2]-dx;
                    x_lright[n+1][2]=0;
                    x_v[n+1][2]=dx-x_lleft[n+1][2]-x_lright[n+1][2];
                    v_lleft[n+1][2]=v_lleft[n][2];
                    v_lright[n+1][2]=0;
                    T_lleft[n+1][2]=T_lleft[n][2];
                    T_lright[n+1][2]=0;
				}
				else//原来为气-液
				{
					sort[n+1]=8;	//现在变成液-气-液
					x_lleft[n+1][2]=x_lleft[n][2]-dx;
					x_v[n+1][2]=dx-x_lleft[n+1][2]-x_lright[n+1][2];
					v_lleft[n+1][2]=v_lleft[n][2];
					T_lleft[n+1][2]=T_lleft[n][2];
				}
				sort[n]=1;	//原右界面变成液
                x_lleft[n][2]=dx;
                x_lright[n][2]=dx;
                x_v[n][2]=0;
                v_lright[n][2]=v_lleft[n][2];
                T_lright[n][2]=T_lleft[n][2];
                ro_v[n][2]=0;
                T_v[n][2]=0;
				T_sat[n][2] = 0;
                P_v[n][2]=0;
			}
			/*右界面向左偏移一个控制体，n、n-1合并考虑*/
			else												
			{
				if (n==1)	//原右界面在第1个控制体，现移动至最后一个控制体
				{
					x_lleft[length][2]=dx+x_lleft[n][2];
                    x_lright[length][2]=0;
                    x_v[length][2]=dx-x_lleft[length][2]-x_lright[length][2];
                    v_lright[length][2]=0;
                    T_lright[length][2]=0;
                    sort[length]=3;	//液-气
				}
				else//原右界面在第2个控制体或者更右
				{
					x_lleft[n-1][2]=dx+x_lleft[n][2];
                    x_lright[n-1][2]=0;
                    x_v[n-1][2]=dx-x_lleft[n-1][2]-x_lright[n-1][2];
                    v_lright[n-1][2]=0;
                    T_lright[n-1][2]=0;
                    sort[n-1]=3;	//液-气
				}
				/*右界面液塞的能量方程*/
				if (n==1)	//原右界面在第1个控制体，现移动至最后一个控制体
				{
					T_lleft[length][2]=1/(x_lleft[length][2]*c_pl_left[length])*(x_lleft[n][1]*T_lleft[n][1]*c_pl_left[n]+x_lleft[length][1]*T_lleft[length][1]*c_pl_left[length]
						+dt/(A_l*ro_l)*(v_lleft[length][1]*A_l*ro_l*c_pl_left[length]*(x_lright[length-1][1]*T_lleft[length][1]+x_lleft[length][1]*T_lright[length-1][1])/(x_lleft[length][1]+x_lright[length-1][1])
						+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]+h_lleft[length]*(T_w[length][1]-T_lleft[length][1])*pi*d*x_lleft[length][1]
						-2*A_l*(T_lleft[length][1]-T_lright[length-1][1])/(x_lleft[length][1]/lamt_left[length]+x_lright[length-1][1]/lamt_right[length-1])
						+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));
					if (T_lleft[length][2]>=T_max)
					{
						T_lleft[length][2]=T_max;
					}
					if (T_lleft[length][2]<=T_min)
					{
						T_lleft[length][2]=T_min;
					}
				}
				else if (n==2)	//原右界面在第2个控制体，现移动至第1个控制体
				{
					T_lleft[n-1][2]=1/(x_lleft[n-1][2]*c_pl_left[n-1])*(x_lleft[n][1]*T_lleft[n][1]*c_pl_left[n]+x_lleft[n-1][1]*T_lleft[n-1][1]*c_pl_left[n-1]
						+dt/(A_l*ro_l)*(v_lleft[n-1][1]*A_l*ro_l*c_pl_left[n-1]*(x_lright[length][1]*T_lleft[n-1][1]+x_lleft[n-1][1]*T_lright[length][1])/(x_lleft[n-1][1]+x_lright[length][1])
						+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]+h_lleft[n-1]*(T_w[n-1][1]-T_lleft[n-1][1])*pi*d*x_lleft[n-1][1]
						-2*A_l*(T_lleft[n-1][1]-T_lright[length][1])/(x_lleft[n-1][1]/lamt_left[n-1]+x_lright[length][1]/lamt_right[length])
						+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));
					if (T_lleft[n-1][2]>=T_max)
					{
						T_lleft[n-1][2]=T_max;
					}
					if (T_lleft[n-1][2]<=T_min)
					{
						T_lleft[n-1][2]=T_min;
					}
				}
				else//原右界面在第3个控制体或者更右，现移动至第2个控制体或更右
				{
					T_lleft[n-1][2]=1/(x_lleft[n-1][2]*c_pl_left[n-1])*(x_lleft[n][1]*T_lleft[n][1]*c_pl_left[n]+x_lleft[n-1][1]*T_lleft[n-1][1]*c_pl_left[n-1]
						+dt/(A_l*ro_l)*(v_lleft[n-1][1]*A_l*ro_l*c_pl_left[n-1]*(x_lright[n-2][1]*T_lleft[n-1][1]+x_lleft[n-1][1]*T_lright[n-2][1])/(x_lleft[n-1][1]+x_lright[n-2][1])
						+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]+h_lleft[n-1]*(T_w[n-1][1]-T_lleft[n-1][1])*pi*d*x_lleft[n-1][1]
						-2*A_l*(T_lleft[n-1][1]-T_lright[n-2][1])/(x_lleft[n-1][1]/lamt_left[n-1]+x_lright[n-2][1]/lamt_right[n-2])
						+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));
					if (T_lleft[n-1][2]>=T_max)
					{
						T_lleft[n-1][2]=T_max;
					}
					if (T_lleft[n-1][2]<=T_min)
					{
						T_lleft[n-1][2]=T_min;
					}
				}
				if (sort[n]==8)	//原右界面为液-气-液
				{
					sort[n]=4;	//现变为气-液
                    x_lleft[n][2]=0;
                    x_v[n][2]=dx-x_lleft[n][2]-x_lright[n][2];
                    T_lleft[n][2]=0;
                    v_lleft[n][2]=0;
				}
				else//原右界面为液-气
				{
					sort[n]=2;	//现变为气
                    T_lleft[n][2]=0;
                    T_lright[n][2]=0;
                    v_lleft[n][2]=0;
                    v_lright[n][2]=0;
                    x_lleft[n][2]=0;
                    x_lright[n][2]=0;
                    x_v[n][2]=dx;
				}
			}

			/*%%%%%%%%%%%%%%%%%%%%%当n_liquid!=1时的讨论%%%%%%%%%%%%%%%%%%%%%%%%%*/
			n_liquid=2;
			while (n_liquid<=n_liquid_total)
			{
				dm_left[n_liquid]=-Tr_v[n_liquid-1]/2*dt;	//左侧气泡分一半质量
				dl_left[n_liquid]=dm_left[n_liquid]/A_l/ro_l;	//左界面移动距离
				dm_right[n_liquid]=-Tr_v[n_liquid]/2*dt;	//右边气泡分一半质量
				dl_right[n_liquid]=dm_right[n_liquid]/A_l/ro_l;    //右界面移动距离
				n=n_trans[2*n_liquid-2];	//左界面
				G_force[n_liquid]=0;
				l[n_liquid]=0;
				C_ll=0;
				K=0;
             
				G_force[n_liquid]=G_force[n_liquid]+ro_l*g*c_gravity[n]*x_lright[n][1];	//各控制体重力累加
				l[n_liquid]=l[n_liquid]+x_lright[n][1];	//各控制体液体长度用x_lright[n][1]表示
				Re_lright[n]=d*ro_l*fabs(v_lright[n][1])/mu_right[n];
				if (Re_lright[n]==0)
					{C_l[n]=0;}
				else if (Re_lright[n]<=1180)
					{C_l[n]=16/Re_lright[n];}
				else
					{C_l[n]=0.078*pow(Re_lright[n],-0.2);}
				if(Re_lright[n]>0&&(n%210>=101)&&(n%210<=210))	//弯头处的阻力损失
					{K=1000/Re_lright[n]+0.1*(1+4/pow(2.3/25.4,0.3));
					C_l[n]=C_l[n]+K/4*2.3/110;}
				if (v_lleft[n][1]>0)
					{C_l[n]=-C_l[n];}
				C_ll=C_ll+C_l[n]*x_lright[n][1];	//阻力系数先加权累加

				n=n+1;
				while (n<=n_trans[2*n_liquid-1])	//从左界面右边控制体到右界面左边控制体
				{
					G_force[n_liquid]=G_force[n_liquid]+ro_l*g*c_gravity[n]*x_lleft[n][1];	//重力累加
					l[n_liquid]=l[n_liquid]+x_lleft[n][1];	//液体长度累加
					Re_lleft[n]=d*ro_l*fabs(v_lleft[n][1])/mu_left[n];
					if (Re_lleft[n]==0)
						{C_l[n]=0;}
					else if (Re_lleft[n]<=1180)
						{C_l[n]=16/Re_lleft[n];}
					else
						{C_l[n]=0.078*pow(Re_lleft[n],-0.2);}
					if(Re_lright[n]>0&&(n%210>=101)&&(n%210<=210))	//弯头处的阻力损失
					{K=1000/Re_lleft[n]+0.1*(1+4/pow(2.3/25.4,0.3));
					C_l[n]=C_l[n]+K/4*2.3/110;}
					if (v_lleft[n][1]>0)
						{C_l[n]=-C_l[n];}
					C_ll=C_ll+C_l[n]*x_lleft[n][1];	//阻力系数加权累加
					n=n+1;
				}
				C_ll=C_ll/l[n_liquid];	//前面加权累加，此处平均

				/*液塞的动量方程，得到右界面速度*/
				v_lleft[n_trans[2*n_liquid-1]][2]=1/(l[n_liquid]+dl_left[n_liquid]+dl_right[n_liquid]-dl_vl_left[n_liquid]-dl_vl_right[n_liquid])
					*(v_lleft[n_trans[2*n_liquid-1]][1]*l[n_liquid]+dt/ro_l*(G_force[n_liquid]+P_v[n_trans[2*n_liquid-2]][1]-P_v[n_trans[2*n_liquid-1]][1]
					+2*C_ll*ro_l/d*l[n_liquid]*pow(v_lleft[n_trans[2*n_liquid-1]][1],2)));

				/*将计算出的v_l值赋给所有同一液塞的控制体*/
				n=n_trans[2*n_liquid-2];	//左界面
				while (n<=n_trans[2*n_liquid-1]-1)	//从左界面到右界面左边
				{
					if (n==n_trans[2*n_liquid-2])	//左界面
					{
						v_lright[n][2]=v_lleft[n_trans[2*n_liquid-1]][2];
					}
					else//左界面以右
					{
						v_lleft[n][2]=v_lleft[n_trans[2*n_liquid-1]][2]; 
						v_lright[n][2]=v_lleft[n_trans[2*n_liquid-1]][2];
					}
					n=n+1;
				}

				/**下面判断是否有位置的偏移（液塞的起点或终点控制体位置改变）**/
				/*左界面*/
				n=n_trans[2*n_liquid-2];
				x_lright[n][2]=x_lright[n][1]-v_lright[n][1]*dt+dl_left[n_liquid]-dl_vl_left[n_liquid];	//左界面新的液体长度
				if (sort[n]!=8)	//左界面为气-液
					{x_lleft[n][2]=x_lleft[n][1];}
				/*若左界面移动超过一个控制体，跳出循环while (n_liquid<=n_liquid_total)，之后还要跳出外层循环*/
				if (x_lright[n][2]<=-dx||x_lright[n][2]>=2*dx)
				{
					codenumber=2;
					break;
				}
				/*左界面无偏移，仍在原控制体*/
				if (x_lright[n][2]>0&&x_lright[n][2]<dx)
				{
					/*%%液体的能量方程%%%*/
					T_lright[n][2]=1/x_lright[n][2]*(x_lright[n][1]*T_lright[n][1]+dt/(c_pl_right[n]*ro_l*A_l)
						*(-v_lright[n][1]*A_l*ro_l*c_pl_right[n]*(x_lleft[n+1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[n+1][1])/(x_lright[n][1]+x_lleft[n+1][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]
						-2*A_l*(T_lright[n][1]-T_lleft[n+1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[n+1][1]/lamt_left[n+1])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));
					if (T_lright[n][2]>=T_max)
					{
						T_lright[n][2]=T_max;
					}
					if (T_lright[n][2]<=T_min)
					{
						T_lright[n][2]=T_min;
					}
					x_v[n][2]=dx-x_lright[n][2]-x_lleft[n][2];
				}
				/*左界面向右移动一个控制体，将n、n+1合并考虑*/
				else if (x_lright[n][2]<=0)								
				{
					x_lright[n+1][2]=dx+x_lright[n][2];
                    x_lleft[n+1][2]=0;
                    x_v[n+1][2]=dx-x_lright[n+1][2]-x_lleft[n+1][2];
                    T_lleft[n+1][2]=0;
                    v_lleft[n+1][2]=0;
                    sort[n+1]=4;	//气-液
                    jump=1;	//左界面右移标记

					/*左界面液体的能量方程*/
                    T_lright[n+1][2]=1/(x_lright[n+1][2]*c_pl_right[n+1])*(x_lright[n][1]*T_lright[n][1]*c_pl_right[n]+x_lright[n+1][1]*T_lright[n+1][1]*c_pl_right[n+1]
						+dt/(A_l*ro_l)*(-v_lright[n+1][1]*A_l*ro_l*c_pl_right[n+1]*(x_lleft[n+2][1]*T_lright[n+1][1]+x_lright[n+1][1]*T_lleft[n+2][1])/(x_lright[n+1][1]+x_lleft[n+2][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]+h_lright[n+1]*(T_w[n+1][1]-T_lright[n+1][1])*pi*d*x_lright[n+1][1]
						-2*A_l*(T_lright[n+1][1]-T_lleft[n+2][1])/(x_lright[n+1][1]/lamt_right[n+1]+x_lleft[n+2][1]/lamt_left[n+2])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));
					if (T_lright[n+1][2]>=T_max)
					{
						T_lright[n+1][2]=T_max;
					}
					if (T_lright[n+1][2]<=T_min)
					{
						T_lright[n+1][2]=T_min;
					}
					if (sort[n]!=8)	//原左界面为气-液
					{
						T_lleft[n][2]=0;
                        T_lright[n][2]=0;
                        v_lleft[n][2]=0;
                        v_lright[n][2]=0;
                        x_lleft[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx;
                        sort[n]=2;	//现变成气
					}
					else//原左界面为液-气-液
					{
						sort[n]=3;	//现变为液-气
                        T_lright[n][2]=0;
                        v_lright[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx-x_lleft[n][2];
					}
				}
				/*左界面向左移动一个控制体，将n、n-1合并考虑*/
				else															
				{
					/*左界面液体的能量方程*/
                    T_lright[n][2]=1/x_lright[n][2]*(x_lright[n][1]*T_lright[n][1]+dt/(c_pl_right[n]*ro_l*A_l)
						*(-v_lright[n][1]*A_l*ro_l*c_pl_right[n]*(x_lleft[n+1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[n+1][1])/(x_lright[n][1]+x_lleft[n+1][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]
						-2*A_l*(T_lright[n][1]-T_lleft[n+1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[n+1][1]/lamt_left[n+1])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));
					if (T_lright[n][2]>=T_max)
					{
						T_lright[n][2]=T_max;
					}
					if (T_lright[n][2]<=T_min)
					{
						T_lright[n][2]=T_min;
					}

					if (sort[n-1]!=3)	//原左界面左侧为气
					{
						x_lright[n-1][2]=x_lright[n][2]-dx; 
                        x_lleft[n-1][2]=0;
                        x_v[n-1][2]=dx-x_lright[n-1][2]-x_lleft[n-1][2];
                        v_lright[n-1][2]=v_lright[n][2];
                        v_lleft[n-1][2]=0;
                        T_lright[n-1][2]=T_lright[n][2];
                        T_lleft[n-1][2]=0;
                        sort[n-1]=4;	//现变为气-液
                          
                        x_lright[n][2]=dx;
                        x_lleft[n][2]=dx;
                        x_v[n][2]=0;
                        v_lleft[n][2]=v_lright[n][2];
                        T_lleft[n][2]=T_lright[n][2];
                        ro_v[n][2]=0;
                        T_v[n][2]=0;
						T_sat[n][2] = 0;
                        P_v[n][2]=0;
                        sort[n]=1;	//原左界面变为液
					}
					else//原左界面为液-气
					{
						x_lright[n-1][2]=x_lright[n][2]-dx; 
                        x_v[n-1][2]=dx-x_lright[n-1][2]-x_lleft[n-1][2];
                        v_lright[n-1][2]=v_lright[n][2];
                        T_lright[n-1][2]=T_lright[n][2];
                        sort[n-1]=8;	//现变为液-气-液
                          
                        x_lright[n][2]=dx;
                        x_lleft[n][2]=dx;
                        x_v[n][2]=0;
                        v_lleft[n][2]=v_lright[n][2];
                        T_lleft[n][2]=T_lright[n][2];
                        ro_v[n][2]=0;
                        T_v[n][2]=0;
						T_sat[n][2] = 0;
						P_v[n][2]=0;
						sort[n]=1;	//原左界面变为液
					}
				}

				/*液塞的中间部分*/
				if (jump==1)//左界面右移
				{
					n=n_trans[2*n_liquid-2]+2;	//定位到现左界面右边控制体
					while (n<=n_trans[2*n_liquid-1]-1)
					{
						T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
							*((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
							-(x_lleft[n+1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[n+1][1])/(x_lleft[n][1]+x_lleft[n+1][1]))
							+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
							-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
							-2*A_l*(T_lleft[n][1]-T_lleft[n+1][1])/(x_lleft[n][1]/lamt_left[n]+x_lleft[n+1][1]/lamt_left[n+1]));
						if (T_lleft[n][2]>=T_max)
						{
							T_lleft[n][2]=T_max;
						}
						if (T_lleft[n][2]<=T_min)
						{
							T_lleft[n][2]=T_min;
						}
						T_lright[n][2]=T_lleft[n][2];
						x_lleft[n][2]=dx;
						x_lright[n][2]=dx;
						x_v[n][2]=0;
						n=n+1;
					}
				}
				else//左界面左移或仍在原控制体
				{
					n=n_trans[2*n_liquid-2]+1;	//定位到原左界面右边控制体
					while (n<=n_trans[2*n_liquid-1]-1)
					{
						T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
							*((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
							-(x_lleft[n+1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lleft[n+1][1])/(x_lleft[n][1]+x_lleft[n+1][1]))
							+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
							-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
							-2*A_l*(T_lleft[n][1]-T_lleft[n+1][1])/(x_lleft[n][1]/lamt_left[n]+x_lleft[n+1][1]/lamt_left[n+1]));
						if (T_lleft[n][2]>=T_max)
						{
							T_lleft[n][2]=T_max;
						}
						if (T_lleft[n][2]<=T_min)
						{
							T_lleft[n][2]=T_min;
						}
						T_lright[n][2]=T_lleft[n][2];
						x_lleft[n][2]=dx;
						x_lright[n][2]=dx;
						x_v[n][2]=0;
						n=n+1;
					}
				}
				jump=0;	//重置左界面移动标记

				/*右界面*/
				n=n_trans[2*n_liquid-1];
				x_lleft[n][2]=x_lleft[n][1]+v_lleft[n][1]*dt+dl_right[n_liquid]-dl_vl_right[n_liquid];	//右界面控制体新的液体长度
				if (sort[n]!=8)	//右界面为液-气
					{x_lright[n][2]=x_lright[n][1];}
				/*若右界面移动超过一个控制体，跳出循环while (n_liquid<=n_liquid_total)，之后还要跳出外层循环*/
				if (x_lleft[n][2]<=-dx||x_lleft[n][2]>=2*dx)
				{
					codenumber=2;
					break;
				}
				/*右界面没有移出原控制体*/
				if (x_lleft[n][2]>0&&x_lleft[n][2]<dx)
				{
					/*右界面液塞的能量方程*/
					T_lleft[n][2]=1/x_lleft[n][2]*(x_lleft[n][1]*T_lleft[n][1]+dt/(c_pl_left[n]*ro_l*A_l)
						*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]*(x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lleft[n][1]+x_lright[n-1][1])
						+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
						-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
						+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));
					if (T_lleft[n][2]>=T_max)
					{
						T_lleft[n][2]=T_max;
					}
					if (T_lleft[n][2]<=T_min)
					{
						T_lleft[n][2]=T_min;
					}
					x_v[n][2]=dx-x_lleft[n][2]-x_lright[n][2];
				}
				/*右界面向右移动一个控制体，n、n+1合并考虑*/
				else if (x_lleft[n][2]>=dx&&x_lleft[n][2]<2*dx)
				{
					/*右界面液塞的能量方程*/
                    T_lleft[n][2]=1/x_lleft[n][2]*(x_lleft[n][1]*T_lleft[n][1]+dt/(c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
						*(x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lleft[n][1]+x_lright[n-1][1])
						+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
						-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
						+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));
                    if (T_lleft[n][2]>=T_max)
					{
						T_lleft[n][2]=T_max;
					}
					if (T_lleft[n][2]<=T_min)
					{
						T_lleft[n][2]=T_min;
					}

					if (sort[n+1]!=4)	//原右界面右边控制体为气
					{
						sort[n+1]=3;	//现变成液-气
                        x_lleft[n+1][2]=x_lleft[n][2]-dx;
                        x_lright[n+1][2]=0;
                        x_v[n+1][2]=dx-x_lleft[n+1][2]-x_lright[n+1][2];
                        v_lleft[n+1][2]=v_lleft[n][2];
                        v_lright[n+1][2]=0;
                        T_lleft[n+1][2]=T_lleft[n][2];
                        T_lright[n+1][2]=0;
					}
					else//原右界面右边控制体为气-液
					{
						sort[n+1]=8;	//现变成液-气-液
                        x_lleft[n+1][2]=x_lleft[n][2]-dx;
                        x_v[n+1][2]=dx-x_lleft[n+1][2]-x_lright[n+1][2];
                        v_lleft[n+1][2]=v_lleft[n][2];
                        T_lleft[n+1][2]=T_lleft[n][2];
					}
					sort[n]=1;	//原右界面变成液
                    x_lleft[n][2]=dx;
                    x_lright[n][2]=dx;
                    x_v[n][2]=0;
                    v_lright[n][2]=v_lleft[n][2];
                    T_lright[n][2]=T_lleft[n][2];
                    ro_v[n][2]=0;
                    T_v[n][2]=0;
					T_sat[n][2] = 0;
                    P_v[n][2]=0;
				}
				/*右界面向左移动一个控制体，n、n-1合并考虑*/
				else													
				{
					x_lleft[n-1][2]=dx+x_lleft[n][2];
                    x_lright[n-1][2]=0;
                    x_v[n-1][2]=dx-x_lleft[n-1][2]-x_lright[n-1][2];
                    v_lright[n-1][2]=0;
                    T_lright[n-1][2]=0;
                    sort[n-1]=3;	//现右界面液-气
					
					/*右界面液体的能量方程*/
                    T_lleft[n-1][2]=1/(x_lleft[n-1][2]*c_pl_left[n-1])*(x_lleft[n][1]*T_lleft[n][1]*c_pl_left[n]+x_lleft[n-1][1]*T_lleft[n-1][1]*c_pl_left[n-1]
						+dt/(A_l*ro_l)*(v_lleft[n-1][1]*A_l*ro_l*c_pl_left[n-1]*(x_lright[n-2][1]*T_lleft[n-1][1]+x_lleft[n-1][1]*T_lright[n-2][1])/(x_lleft[n-1][1]+x_lright[n-2][1])
						+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]+h_lleft[n-1]*(T_w[n-1][1]-T_lleft[n-1][1])*pi*d*x_lleft[n-1][1]
						-2*A_l*(T_lleft[n-1][1]-T_lright[n-2][1])/(x_lleft[n-1][1]/lamt_left[n-1]+x_lright[n-2][1]/lamt_right[n-2])
						+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));
					if (T_lleft[n-1][2]>=T_max)
					{
						T_lleft[n-1][2]=T_max;
					}
					if (T_lleft[n-1][2]<=T_min)
					{
						T_lleft[n-1][2]=T_min;
					}

					if (sort[n]==8)	//原右界面为液-气-液
					{
						sort[n]=4;	//现变为气-液
                        x_lleft[n][2]=0;
                        x_v[n][2]=dx-x_lleft[n][2]-x_lright[n][2];
                        T_lleft[n][2]=0;
                        v_lleft[n][2]=0;
					}
					else//原右界面为液-气
					{
						sort[n]=2;	//现变为气
                        T_lleft[n][2]=0;
                        T_lright[n][2]=0;
                        v_lleft[n][2]=0;
                        v_lright[n][2]=0;
                        x_lleft[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx;
					}
				}

				n_liquid=n_liquid+1;
			}

			/*若左、右界面超出范围，跳出外层循环*/
			if (codenumber==2)
				{break;}
		}

		/******对于起始点状态的分类讨论：第三部分：发生异常，出现气-液-气之类*********/
		else
		{
			fid0=fopen("state.txt","a+");	//按照a+模式打开一个名叫state.txt的文件，返回state.txt的文件指针给fid0,异常则返回NULL
			fprintf(fid0,"Strange,sort[1] impossible, n=%d, i=%d",n,i);	//将“Strange,sort[1] impossible”以及n和i写入由fid0指出的文件
			fclose(fid0);	//关闭由fid0指出的文件,返回操作结果，0或EOF
		}

		/*%%%%%%%%%%%%%%%%%%%有气泡通过x=0导致的气泡序号的变化，而液塞序号不变%%%%%%%%%%%%%%%%%%%%%%%%*/
		double P_vb_tem[100][3];	//气泡序号变化后，用来放新值的容器
		double T_vb_tem[100][3];	//气泡序号变化后，用来放新值的容器
		double T_sat_vb_tem[100][3];	//气泡序号变化后，用来放新值的容器
		double ro_vb_tem[100][3];	//气泡序号变化后，用来放新值的容器
		double V_bubble_tem[100][3];	//气泡序号变化后，用来放新值的容器
		double m_bubble_tem[100][3];	//气泡序号变化后，用来放新值的容器
		double ro_vb[100][3];	//整个气泡密度

		if (cross_negative==1)	//气泡序号减1
		{
			n_bubble=1;
			while (n_bubble<=n_bubble_total-1)
			{
				P_vb_tem[n_bubble][1]=P_vb[n_bubble+1][1];
				T_vb_tem[n_bubble][1]=T_vb[n_bubble+1][1];
				T_sat_vb_tem[n_bubble][1] = T_sat_vb[n_bubble + 1][1];
				ro_vb_tem[n_bubble][1]=ro_vb[n_bubble+1][1];
				V_bubble_tem[n_bubble][1]=V_bubble[n_bubble+1][1];
				m_bubble_tem[n_bubble][1]=m_bubble[n_bubble+1][1];
				n_bubble=n_bubble+1;
			}
			P_vb_tem[n_bubble_total][1]=P_vb[1][1];	//第1个气泡，变成最后一个
			T_vb_tem[n_bubble_total][1]=T_vb[1][1];
			T_sat_vb_tem[n_bubble_total][1]=T_sat_vb[1][1];
			ro_vb_tem[n_bubble_total][1]=ro_vb[1][1];
			V_bubble_tem[n_bubble_total][1]=V_bubble[1][1];
			m_bubble_tem[n_bubble_total][1]=m_bubble[1][1];

			n_bubble=1;
			while (n_bubble<=n_bubble_total)	//序号变化后的状态参数
			{
				P_vb[n_bubble][1]=P_vb_tem[n_bubble][1];
				T_vb[n_bubble][1]=T_vb_tem[n_bubble][1];
				T_sat_vb[n_bubble][1] = T_sat_vb_tem[n_bubble][1];
				ro_vb[n_bubble][1]=ro_vb_tem[n_bubble][1];
				V_bubble[n_bubble][1]=V_bubble_tem[n_bubble][1];
				m_bubble[n_bubble][1]=m_bubble_tem[n_bubble][1];
				n_bubble=n_bubble+1;
			}
		}

		if (cross_positive==1)	//气泡序号加1
		{
			n_bubble=1;
			while (n_bubble<=n_bubble_total-1)
			{
				P_vb_tem[n_bubble+1][1]=P_vb[n_bubble][1];
				T_vb_tem[n_bubble+1][1]=T_vb[n_bubble][1];
				T_sat_vb_tem[n_bubble + 1][1] = T_sat_vb[n_bubble][1];
				ro_vb_tem[n_bubble+1][1]=ro_vb[n_bubble][1];
				V_bubble_tem[n_bubble+1][1]=V_bubble[n_bubble][1];
				m_bubble_tem[n_bubble+1][1]=m_bubble[n_bubble][1];
				n_bubble=n_bubble+1;
			}
			P_vb_tem[1][1]=P_vb[n_bubble_total][1];	//最后一个气泡变成第1个气泡
			T_vb_tem[1][1]=T_vb[n_bubble_total][1];
			T_sat_vb_tem[1][1] = T_sat_vb[n_bubble_total][1];
			ro_vb_tem[1][1]=ro_vb[n_bubble_total][1];
			V_bubble_tem[1][1]=V_bubble[n_bubble_total][1];
			m_bubble_tem[1][1]=m_bubble[n_bubble_total][1];

			n_bubble=1;
			while (n_bubble<=n_bubble_total)
			{
				P_vb[n_bubble][1]=P_vb_tem[n_bubble][1];
				T_vb[n_bubble][1]=T_vb_tem[n_bubble][1];
				T_sat_vb[n_bubble][1] = T_sat_vb_tem[n_bubble][1];
				ro_vb[n_bubble][1]=ro_vb_tem[n_bubble][1];
				V_bubble[n_bubble][1]=V_bubble_tem[n_bubble][1];
				m_bubble[n_bubble][1]=m_bubble_tem[n_bubble][1];
				n_bubble=n_bubble+1;
			}		
		}


		/*%%%%对出现sort(n)==9(气泡左右界面缩到同一控制体：通常导致1113111或1114111)%%%*/

		int n_sort9;	//sort[n]==9的控制体个数
		int n_bubble_tem;
		n=1;
		if ((sort[1]==3||sort[1]==4)&&sort[length]*sort[2]==1)	//气泡左右界面都在第1个控制体
			{sort[1]=9;}
		n=length;
		if ((sort[length]==3||sort[length]==4)&&sort[length-1]*sort[1]==1)	//气泡左右界面都在最后一个控制体
			{sort[length]=9;}
		n=2;
		while (n<=length-1)
		{
			if ((sort[n]==3||sort[n]==4)&&sort[n-1]*sort[n+1]==1)
				{sort[n]=9;}
			n=n+1;
		}
        
		/*若有两个及以上sort[n]==9的点则特殊输出*/
		n_sort9=0;
		n=1;
		while (n<=length)
		{
			if (sort[n]==9)
				{n_sort9=n_sort9+1;}
			n=n+1;
		}
		if (n_sort9>=2)
		{
			fid0=fopen("state.txt","a+");	//按照a+模式打开一个名叫state.txt的文件，返回state.txt的文件指针给fid0,异常则返回NULL
			fprintf(fid0,"Special, Codenumber=%d, n=%d, i=%d, n_sort9=%d",codenumber,n,i,n_sort9);	//将“”写入由fid0指出的文件，此时codenumber==？？？
			fclose(fid0);	//关闭由fid0指出的文件,返回操作结果，0或EOF
		}
		codenumber=0;

		/*消灭sort[n]==9，气泡序号正式变化*/
		if (n_sort9!=0)
		{
			if (sort[1]==1)	//第1个控制体为液
			{
				n_bubble=1;	//初始为1，给sort[n]==9一个气泡
				n=2;
				transfer_point_number=0;
				while (n<=length)
				{
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=transfer_point_number+1;
						if ((transfer_point_number%2)==0)	//整除2，余数0
							{n_bubble=n_bubble+1;}
					}
					if (sort[n]==8)
					{
						transfer_point_number=transfer_point_number+2;
						n_bubble=n_bubble+1;
					}
					if (sort[n]==9)	//把sort[n]=9改造成sort[n]=1，气泡序号减小1
					{
						x_lleft[n][2]=dx;
						x_lright[n][2]=dx;
						x_v[n][2]=0;
						if (n==length)
						{
							v_lleft[n][2]=0.5*(v_lleft[1][2]+v_lright[n-1][2]);
							v_lright[n][2]=v_lleft[n][2];
							T_lleft[n][2]=0.5*(T_lleft[1][2]+T_lright[n-1][2]);
							T_lright[n][2]=T_lleft[n][2];
						}
						else
						{
							v_lleft[n][2]=0.5*(v_lleft[n+1][2]+v_lright[n-1][2]);
							v_lright[n][2]=v_lleft[n][2];
							T_lleft[n][2]=0.5*(T_lleft[n+1][2]+T_lright[n-1][2]);
							T_lright[n][2]=T_lleft[n][2];
						}
						P_v[n][2]=0;
						T_v[n][2]=0;
						T_sat[n][2] = 0;
						ro_v[n][2]=0;
						sort[n]=1;

						n_bubble_tem=n_bubble;
						while (n_bubble_tem<=n_bubble_total-1)
						{
							P_vb[n_bubble_tem][1]=P_vb[n_bubble_tem+1][1];
							T_vb[n_bubble_tem][1]=T_vb[n_bubble_tem+1][1];
							T_sat_vb[n_bubble_tem][1] = T_sat_vb[n_bubble_tem + 1][1];
							ro_vb[n_bubble_tem][1]=ro_vb[n_bubble_tem+1][1];
							V_bubble[n_bubble_tem][1]=V_bubble[n_bubble_tem+1][1];
							m_bubble[n_bubble_tem][1]=m_bubble[n_bubble_tem+1][1];
							n_bubble_tem=n_bubble_tem+1;
						}
					}
					n=n+1;
				}
			}
			if (sort[1]==2)	//第1个控制体为气
			{
				n_bubble=1;	//提前给sort[n]==9一个气泡数量
				n=2;
				transfer_point_number=0;
				while (n<=length)
				{
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=transfer_point_number+1;
						if ((transfer_point_number%2)!=0)	//不能整除2
							{n_bubble=n_bubble+1;}
					}
					if (sort[n]==8)
					{
						transfer_point_number=transfer_point_number+2;
						n_bubble=n_bubble+1;
					}
					if (sort[n]==9)	//改造成sort[n]==1，气泡序号减小1
					{
						x_lleft[n][2]=dx;
						x_lright[n][2]=dx;
						x_v[n][2]=0;
						if (n==length)
						{
							v_lleft[n][2]=0.5*(v_lleft[1][2]+v_lright[n-1][2]);
							v_lright[n][2]=v_lleft[n][2];
							T_lleft[n][2]=0.5*(T_lleft[1][2]+T_lright[n-1][2]);
							T_lright[n][2]=T_lleft[n][2];
						}
						else
						{
							v_lleft[n][2]=0.5*(v_lleft[n+1][2]+v_lright[n-1][2]);
							v_lright[n][2]=v_lleft[n][2];
							T_lleft[n][2]=0.5*(T_lleft[n+1][2]+T_lright[n-1][2]);
							T_lright[n][2]=T_lleft[n][2];
						}
						P_v[n][2]=0;
						T_v[n][2]=0;
						T_sat[n][2]=0;
						ro_v[n][2]=0;
						sort[n]=1;

						n_bubble_tem=n_bubble;
						while (n_bubble_tem<=n_bubble_total-1)	//sort[n]==9之后的气泡序号减小1
						{
							P_vb[n_bubble_tem][1]=P_vb[n_bubble_tem+1][1];
							T_vb[n_bubble_tem][1]=T_vb[n_bubble_tem+1][1];
							T_sat_vb[n_bubble_tem][1]=T_sat_vb[n_bubble_tem+1][1];
							ro_vb[n_bubble_tem][1]=ro_vb[n_bubble_tem+1][1];
							V_bubble[n_bubble_tem][1]=V_bubble[n_bubble_tem+1][1];
							m_bubble[n_bubble_tem][1]=m_bubble[n_bubble_tem+1][1];
							n_bubble_tem=n_bubble_tem+1;
						}
					}
					n=n+1;
				}
			}
			if (sort[1]==3)	//第1个控制体为液-气
			{
				n_bubble=1;
				n=2;
				transfer_point_number=1;
				while (n<=length)
				{
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=transfer_point_number+1;
						if ((transfer_point_number%2)==0)	//整除2
							{n_bubble=n_bubble+1;}
					}
					if (sort[n]==8)
					{
						transfer_point_number=transfer_point_number+2;
						n_bubble=n_bubble+1;
					}
					if (sort[n]==9)//改造成sort[n]==1，气泡序号减小1
					{
						x_lleft[n][2]=dx;
						x_lright[n][2]=dx;
						x_v[n][2]=0;
						if (n==length)
						{
							v_lleft[n][2]=0.5*(v_lleft[1][2]+v_lright[n-1][2]);
							v_lright[n][2]=v_lleft[n][2];
							T_lleft[n][2]=0.5*(T_lleft[1][2]+T_lright[n-1][2]);
							T_lright[n][2]=T_lleft[n][2];
						}
						else
						{
							v_lleft[n][2]=0.5*(v_lleft[n+1][2]+v_lright[n-1][2]);
							v_lright[n][2]=v_lleft[n][2];
							T_lleft[n][2]=0.5*(T_lleft[n+1][2]+T_lright[n-1][2]);
							T_lright[n][2]=T_lleft[n][2];
						}
						P_v[n][2]=0;
						T_v[n][2]=0;
						T_sat[n][2] = 0;
						ro_v[n][2]=0;
						sort[n]=1;

						n_bubble_tem=n_bubble;
						while (n_bubble_tem<=n_bubble_total-1)
						{
							P_vb[n_bubble_tem][1]=P_vb[n_bubble_tem+1][1];
							T_vb[n_bubble_tem][1]=T_vb[n_bubble_tem+1][1];
							T_sat_vb[n_bubble_tem][1] = T_sat_vb[n_bubble_tem + 1][1];
							ro_vb[n_bubble_tem][1]=ro_vb[n_bubble_tem+1][1];
							V_bubble[n_bubble_tem][1]=V_bubble[n_bubble_tem+1][1];
							m_bubble[n_bubble_tem][1]=m_bubble[n_bubble_tem+1][1];
							n_bubble_tem=n_bubble_tem+1;
						}
					}
					n=n+1;
				}
			}
			if (sort[1]==4)	//第1个控制体为气-液
			{
				n_bubble=2;
				n=2;
				transfer_point_number=1;
				while (n<=length)
				{
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=transfer_point_number+1;
						if ((transfer_point_number%2)!=0)	//不能整除2，即为奇数
							{n_bubble=n_bubble+1;}
					}
					if (sort[n]==8)
					{
						transfer_point_number=transfer_point_number+2;
						n_bubble=n_bubble+1;
					}
					if (sort[n]==9)	//改造成sort[n]==1，气泡序号减小1
					{
						x_lleft[n][2]=dx;
						x_lright[n][2]=dx;
						x_v[n][2]=0;
						if (n==length)
						{
							v_lleft[n][2]=0.5*(v_lleft[1][2]+v_lright[n-1][2]);
							v_lright[n][2]=v_lleft[n][2];
							T_lleft[n][2]=0.5*(T_lleft[1][2]+T_lright[n-1][2]);
							T_lright[n][2]=T_lleft[n][2];
						}
						else
						{
							v_lleft[n][2]=0.5*(v_lleft[n+1][2]+v_lright[n-1][2]);
							v_lright[n][2]=v_lleft[n][2];
							T_lleft[n][2]=0.5*(T_lleft[n+1][2]+T_lright[n-1][2]);
							T_lright[n][2]=T_lleft[n][2];
						}
						P_v[n][2]=0;
						T_v[n][2]=0;
						T_sat[n][2] = 0;
						ro_v[n][2]=0;
						sort[n]=1;

						n_bubble_tem=n_bubble;
						while (n_bubble_tem<=n_bubble_total-1)
						{
							P_vb[n_bubble_tem][1]=P_vb[n_bubble_tem+1][1];
							T_vb[n_bubble_tem][1]=T_vb[n_bubble_tem+1][1];
							T_sat_vb[n_bubble_tem][1] = T_sat_vb[n_bubble_tem + 1][1];
							ro_vb[n_bubble_tem][1]=ro_vb[n_bubble_tem+1][1];
							V_bubble[n_bubble_tem][1]=V_bubble[n_bubble_tem+1][1];
							m_bubble[n_bubble_tem][1]=m_bubble[n_bubble_tem+1][1];
							n_bubble_tem=n_bubble_tem+1;
						}
					}
					n=n+1;
				}
			}
			if (sort[1]==9)	//第1个控制体包含一整个气泡
			{
				n_bubble=1;
				n=1;
				transfer_point_number=2;
				x_lleft[n][2]=dx;
				x_lright[n][2]=dx;
				x_v[n][2]=0;
				v_lleft[n][2]=0.5*(v_lleft[n+1][2]+v_lright[length][2]);
				v_lright[n][2]=v_lleft[n][2];
				T_lleft[n][2]=0.5*(T_lleft[n+1][2]+T_lright[length][2]);
				T_lright[n][2]=T_lleft[n][2];
				P_v[n][2]=0;
				T_v[n][2]=0;
				T_sat[n][2] = 0;
				ro_v[n][2]=0;
				sort[n]=1;	//改造成sort[n]==1，气泡序号减小1

				n_bubble_tem=n_bubble;
				while (n_bubble_tem<=n_bubble_total-1)
				{
					P_vb[n_bubble_tem][1]=P_vb[n_bubble_tem+1][1];
					T_vb[n_bubble_tem][1]=T_vb[n_bubble_tem+1][1];
					T_sat_vb[n_bubble_tem][1] = T_sat_vb[n_bubble_tem + 1][1];
					ro_vb[n_bubble_tem][1]=ro_vb[n_bubble_tem+1][1];
					V_bubble[n_bubble_tem][1]=V_bubble[n_bubble_tem+1][1];
					m_bubble[n_bubble_tem][1]=m_bubble[n_bubble_tem+1][1];
					n_bubble_tem=n_bubble_tem+1;
				}

				n=2;
				transfer_point_number=2;	//多余
				n_bubble=2;	//多余，为何不是1？？？
				while (n<=length)
				{
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=transfer_point_number+1;
						if ((transfer_point_number%2)==0)
							{n_bubble=n_bubble+1;}
					}
					if (sort[n]==8)
					{
						transfer_point_number=transfer_point_number+2;
						n_bubble=n_bubble+1;
					}
					if (sort[n]==9)	//改造成sort[n]==1，气泡序号减小1
					{
						x_lleft[n][2]=dx;
						x_lright[n][2]=dx;
						x_v[n][2]=0;
						if (n==length)
						{
							v_lleft[n][2]=0.5*(v_lleft[1][2]+v_lright[n-1][2]);
							v_lright[n][2]=v_lleft[n][2];
							T_lleft[n][2]=0.5*(T_lleft[1][2]+T_lright[n-1][2]);
							T_lright[n][2]=T_lleft[n][2];
						}
						else
						{
							v_lleft[n][2]=0.5*(v_lleft[n+1][2]+v_lright[n-1][2]);
							v_lright[n][2]=v_lleft[n][2];
							T_lleft[n][2]=0.5*(T_lleft[n+1][2]+T_lright[n-1][2]);
							T_lright[n][2]=T_lleft[n][2];
						}
						P_v[n][2]=0;
						T_v[n][2]=0;
						T_sat[n][2] = 0;
						ro_v[n][2]=0;
						sort[n]=1;

						n_bubble_tem=n_bubble;
						while (n_bubble_tem<=n_bubble_total-1)
						{
							P_vb[n_bubble_tem][1]=P_vb[n_bubble_tem+1][1];
							T_vb[n_bubble_tem][1]=T_vb[n_bubble_tem+1][1];
							T_sat_vb[n_bubble_tem][1] = T_sat_vb[n_bubble_tem + 1][1];
							ro_vb[n_bubble_tem][1]=ro_vb[n_bubble_tem+1][1];
							V_bubble[n_bubble_tem][1]=V_bubble[n_bubble_tem+1][1];
							m_bubble[n_bubble_tem][1]=m_bubble[n_bubble_tem+1][1];
							n_bubble_tem=n_bubble_tem+1;
						}
					}
					n=n+1;
				}
			}
		}
		n_sort9=0;

		/*%%%%%%液塞合并位置的确认，新分布下气泡值的重新定义%%%%%%%*/
		if (sort[1]==1)	//第1个控制体为液
		{
			n_bubble=1;
			n=2;
			transfer_point_number=0;
			while (n<=length)
			{
				if (sort[n]==3||sort[n]==4)
				{
					transfer_point_number=transfer_point_number+1;
					if ((transfer_point_number%2)==0)
						{n_bubble=n_bubble+1;}
				}
				if (sort[n]==8)	//液塞在此可能合并
				{
					x_v[n][2]=dx-x_lleft[n][2]-x_lright[n][2];
					if (x_v[n][2]<l_disappear)	//液塞合并
					{
						x_lleft[n][2]=dx;
						x_lright[n][2]=dx;
						x_v[n][2]=0;
						v_lleft[n][2]=0.5*(v_lleft[n][2]+v_lright[n][2]);
						v_lright[n][2]=v_lleft[n][2];
						T_lleft[n][2]=0.5*(T_lleft[n][2]+T_lright[n][2]);
						T_lright[n][2]=T_lleft[n][2];
						P_v[n][2]=0;
						T_v[n][2]=0;
						T_sat[n][2] = 0;
						ro_v[n][2]=0;
						sort[n]=1;	//合并后变成1，此后的气泡序号减小1
                   
						n_bubble_tem=n_bubble;
						while (n_bubble_tem<=n_bubble_total-1)
						{
							P_vb[n_bubble_tem][1]=P_vb[n_bubble_tem+1][1];
							T_vb[n_bubble_tem][1]=T_vb[n_bubble_tem+1][1];
							T_sat_vb[n_bubble_tem][1] = T_sat_vb[n_bubble_tem + 1][1];
							ro_vb[n_bubble_tem][1]=ro_vb[n_bubble_tem+1][1];
							V_bubble[n_bubble_tem][1]=V_bubble[n_bubble_tem+1][1];
							m_bubble[n_bubble_tem][1]=m_bubble[n_bubble_tem+1][1];
							n_bubble_tem=n_bubble_tem+1;
						}
					}
					else  //液塞不合并
					{
						transfer_point_number=transfer_point_number+2;
						if ((transfer_point_number%2)==0)
							{n_bubble=n_bubble+1;}
					}
				}
				n=n+1;
			}
		}
		if (sort[1]==2)	//第1个控制体为气
		{
			n_bubble=1;
			n=2;
			transfer_point_number=0;
			while (n<=length)
			{
				if (sort[n]==3||sort[n]==4)
				{
					transfer_point_number=transfer_point_number+1;
					if ((transfer_point_number%2)!=0)
						{n_bubble=n_bubble+1;}
				}
				if (sort[n]==8)
				{
					x_v[n][2]=dx-x_lleft[n][2]-x_lright[n][2];
					if (x_v[n][2]<l_disappear)
					{
						x_lleft[n][2]=dx;
						x_lright[n][2]=dx;
						x_v[n][2]=0;
						v_lleft[n][2]=0.5*(v_lleft[n][2]+v_lright[n][2]);
						v_lright[n][2]=v_lleft[n][2];
						T_lleft[n][2]=0.5*(T_lleft[n][2]+T_lright[n][2]);
						T_lright[n][2]=T_lleft[n][2];
						P_v[n][2]=0;
						T_v[n][2]=0;
						T_sat[n][2] = 0;
						ro_v[n][2]=0;
						sort[n]=1;	//液塞合并，此后气泡序号减小1
                   
						n_bubble_tem=n_bubble;
						while (n_bubble_tem<=n_bubble_total-1)
						{
							P_vb[n_bubble_tem][1]=P_vb[n_bubble_tem+1][1];
							T_vb[n_bubble_tem][1]=T_vb[n_bubble_tem+1][1];
							T_sat_vb[n_bubble_tem][1] = T_sat_vb[n_bubble_tem + 1][1];
							ro_vb[n_bubble_tem][1]=ro_vb[n_bubble_tem+1][1];
							V_bubble[n_bubble_tem][1]=V_bubble[n_bubble_tem+1][1];
							m_bubble[n_bubble_tem][1]=m_bubble[n_bubble_tem+1][1];
							n_bubble_tem=n_bubble_tem+1;
						}
					}
					else  //液塞不合并
					{
						transfer_point_number=transfer_point_number+2;
						if ((transfer_point_number%2)!=0)
							{n_bubble=n_bubble+1;}
					}
				}
				n=n+1;
			}
		}
		if (sort[1]==3)	//第1个控制体为液-气
		{
			n_bubble=1;
			n=2;
			transfer_point_number=1;
			while (n<=length)
			{
				if (sort[n]==3||sort[n]==4)
				{
					transfer_point_number=transfer_point_number+1;
					if ((transfer_point_number%2)==0)
						{n_bubble=n_bubble+1;}
				}
				if (sort[n]==8)
				{
					x_v[n][2]=dx-x_lleft[n][2]-x_lright[n][2];
					if (x_v[n][2]<l_disappear)	//液塞合并
					{
						x_lleft[n][2]=dx;
						x_lright[n][2]=dx;
						x_v[n][2]=0;
						v_lleft[n][2]=0.5*(v_lleft[n][2]+v_lright[n][2]);
						v_lright[n][2]=v_lleft[n][2];
						T_lleft[n][2]=0.5*(T_lleft[n][2]+T_lright[n][2]);
						T_lright[n][2]=T_lleft[n][2];
						P_v[n][2]=0;
						T_v[n][2]=0;
						T_sat[n][2] = 0;
						ro_v[n][2]=0;
						sort[n]=1;
                   
						n_bubble_tem=n_bubble;
						while (n_bubble_tem<=n_bubble_total-1)
						{
							P_vb[n_bubble_tem][1]=P_vb[n_bubble_tem+1][1];
							T_vb[n_bubble_tem][1]=T_vb[n_bubble_tem+1][1];
							T_sat_vb[n_bubble_tem][1] = T_sat_vb[n_bubble_tem + 1][1];
							ro_vb[n_bubble_tem][1]=ro_vb[n_bubble_tem+1][1];
							V_bubble[n_bubble_tem][1]=V_bubble[n_bubble_tem+1][1];
							m_bubble[n_bubble_tem][1]=m_bubble[n_bubble_tem+1][1];
							n_bubble_tem=n_bubble_tem+1;
						}
					}
					else  //液塞不合并
					{
						transfer_point_number=transfer_point_number+2;
						if ((transfer_point_number%2)==0)
							{n_bubble=n_bubble+1;}
					}
				}
				n=n+1;
			}

		}
		if (sort[1]==4)	//第1个控制体为气-液
		{
			n_bubble=2;
			n=2;
			transfer_point_number=1;
			while (n<=length)
			{
				if (sort[n]==3||sort[n]==4)
				{
					transfer_point_number=transfer_point_number+1;
					if ((transfer_point_number%2)!=0)
						{n_bubble=n_bubble+1;}
				}
				if (sort[n]==8)
				{
					x_v[n][2]=dx-x_lleft[n][2]-x_lright[n][2];
					if (x_v[n][2]<l_disappear)	//液塞合并
					{
						x_lleft[n][2]=dx;
						x_lright[n][2]=dx;
						x_v[n][2]=0;
						v_lleft[n][2]=0.5*(v_lleft[n][2]+v_lright[n][2]);
						v_lright[n][2]=v_lleft[n][2];
						T_lleft[n][2]=0.5*(T_lleft[n][2]+T_lright[n][2]);
						T_lright[n][2]=T_lleft[n][2];
						P_v[n][2]=0;
						T_v[n][2]=0;
						T_sat[n][2] = 0;
						ro_v[n][2]=0;
						sort[n]=1;
                   
						n_bubble_tem=n_bubble;
						while (n_bubble_tem<=n_bubble_total-1)
						{
							P_vb[n_bubble_tem][1]=P_vb[n_bubble_tem+1][1];
							T_vb[n_bubble_tem][1]=T_vb[n_bubble_tem+1][1];
							T_sat_vb[n_bubble_tem][1] = T_sat_vb[n_bubble_tem + 1][1];
							ro_vb[n_bubble_tem][1]=ro_vb[n_bubble_tem+1][1];
							V_bubble[n_bubble_tem][1]=V_bubble[n_bubble_tem+1][1];
							m_bubble[n_bubble_tem][1]=m_bubble[n_bubble_tem+1][1];
							n_bubble_tem=n_bubble_tem+1;
						}
					}
					else  //液塞不合并
					{
						transfer_point_number=transfer_point_number+2;
						if ((transfer_point_number%2)!=0)
							{n_bubble=n_bubble+1;}
					}
				}
				n=n+1;
			}

		}
		if (sort[1]==8)	//第1个控制体为液-气-液
		{
			n_bubble=1;
			x_v[1][2]=dx-x_lleft[1][2]-x_lright[1][2];
			if (x_v[1][2]<l_disappear)	//液塞合并
			{
				x_lleft[1][2]=dx;
				x_lright[1][2]=dx;
				x_v[1][2]=0;
				v_lleft[1][2]=0.5*(v_lleft[1][2]+v_lright[1][2]);
				v_lright[1][2]=v_lleft[1][2];
				T_lleft[1][2]=0.5*(T_lleft[1][2]+T_lright[1][2]);
				T_lright[1][2]=T_lleft[1][2];
				P_v[1][2]=0;
				T_v[1][2]=0;
				T_sat[1][2] = 0;
				ro_v[1][2]=0;
				sort[1]=1;
                   
				n_bubble_tem=n_bubble;
				while (n_bubble_tem<=n_bubble_total-1)
				{
					P_vb[n_bubble_tem][1]=P_vb[n_bubble_tem+1][1];
					T_vb[n_bubble_tem][1]=T_vb[n_bubble_tem+1][1];
					T_sat_vb[n_bubble_tem][1] = T_sat_vb[n_bubble_tem + 1][1];
					ro_vb[n_bubble_tem][1]=ro_vb[n_bubble_tem+1][1];
					V_bubble[n_bubble_tem][1]=V_bubble[n_bubble_tem+1][1];
					m_bubble[n_bubble_tem][1]=m_bubble[n_bubble_tem+1][1];
					n_bubble_tem=n_bubble_tem+1;
				}
			}
			else  //液塞不合并
			{
				transfer_point_number=2;
				n_bubble=n_bubble+1;
			}

			n=2;
			while (n<=length)
			{
				if (sort[n]==3||sort[n]==4)
				{
					transfer_point_number=transfer_point_number+1;
					if ((transfer_point_number%2)==0)
						{n_bubble=n_bubble+1;}
				}
				if (sort[n]==8)
				{
					x_v[n][2]=dx-x_lleft[n][2]-x_lright[n][2];
					if (x_v[n][2]<l_disappear)	//液塞合并
					{
						x_lleft[n][2]=dx;
						x_lright[n][2]=dx;
						x_v[n][2]=0;
						v_lleft[n][2]=0.5*(v_lleft[n][2]+v_lright[n][2]);
						v_lright[n][2]=v_lleft[n][2];
						T_lleft[n][2]=0.5*(T_lleft[n][2]+T_lright[n][2]);
						T_lright[n][2]=T_lleft[n][2];
						P_v[n][2]=0;
						T_v[n][2]=0;
						T_sat[n][2] = 0;
						ro_v[n][2]=0;
						sort[n]=1;
                   
						n_bubble_tem=n_bubble;
						while (n_bubble_tem<=n_bubble_total-1)
						{
							P_vb[n_bubble_tem][1]=P_vb[n_bubble_tem+1][1];
							T_vb[n_bubble_tem][1]=T_vb[n_bubble_tem+1][1];
							T_sat_vb[n_bubble_tem][1] = T_sat_vb[n_bubble_tem + 1][1];
							ro_vb[n_bubble_tem][1]=ro_vb[n_bubble_tem+1][1];
							V_bubble[n_bubble_tem][1]=V_bubble[n_bubble_tem+1][1];
							m_bubble[n_bubble_tem][1]=m_bubble[n_bubble_tem+1][1];
							n_bubble_tem=n_bubble_tem+1;
						}
					}
					else  //液塞不合并
					{
						transfer_point_number=transfer_point_number+2;
						if ((transfer_point_number%2)==0)
							{n_bubble=n_bubble+1;}
					}
				}
				n=n+1;
			}
		}
		/*重新确定气液界面位置和数量*/
		n=1;
		if (sort[n]==3||sort[n]==4)
		{
			transfer_point_number=1;
			n_trans[1]=1;
		}
		else if (sort[n]==8)
		{
			transfer_point_number=2;
            n_trans[1]=1;
            n_trans[2]=1;
		}
		else
		{
			transfer_point_number=0;
		}
		n=2;
		while (n<=length)
		{
			if (sort[n]==3||sort[n]==4)
			{
				transfer_point_number=transfer_point_number+1;
				n_trans[transfer_point_number]=n;                  //使用n_trans[]来记录每个交界面的位置
			}
			if (sort[n]==8)
			{
				transfer_point_number=transfer_point_number+1;
				n_trans[transfer_point_number]=n;                  
				transfer_point_number=transfer_point_number+1;
				n_trans[transfer_point_number]=n;                  //n_trans{]来记录每个交界面的位置
			}
			n=n+1;
		}

		n_trans_total=transfer_point_number;	//记录气液界面总数
		n_bubble_total=n_trans_total/2;	//气泡的数量
		n_liquid_total=n_trans_total/2;  //液塞的数量

		/*%%%%%%%%计算液塞合并之后的液塞速度%%%%%%%%*/
		double momentum;	//单位动量，m2/s
		double l_liquid;	//液塞长度，m

		if (sort[1]==2||sort[1]==4)	//第1个控制体为气或者气-液，此时液塞序号领先气泡序号
		{
			n_liquid=1;
			while (n_liquid<=n_liquid_total)
			{
				n=n_trans[2*n_liquid-1];	//液塞左界面
				momentum=0;
				l_liquid=0;
				momentum=momentum+v_lright[n][2]*x_lright[n][2];	//两边同除以密度和截面积
				l_liquid=l_liquid+x_lright[n][2];	
				n=n+1;
				while (n<=n_trans[2*n_liquid])	//从左界面右侧控制体到右界面
				{
					momentum=momentum+v_lleft[n][2]*x_lleft[n][2];	//动量加权累加
					l_liquid=l_liquid+x_lleft[n][2];
					n=n+1; 
				}
				n=n_trans[2*n_liquid-1];	//左界面
				v_lright[n][2]=momentum/l_liquid;	//新的速度
				n=n+1;
				while (n<=n_trans[2*n_liquid]-1)	//把左界面的速度赋给同一液塞的所有控制体
				{
					v_lright[n][2]=momentum/l_liquid;
					v_lleft[n][2]=momentum/l_liquid;
					n=n+1;
				}
				n=n_trans[2*n_liquid];	//右界面
				v_lleft[n][2]=momentum/l_liquid;

				n_liquid=n_liquid+1;
			}
		}
		if (sort[1]==1||sort[1]==3||sort[1]==8)	//第1个控制体为液、液-气、液-气-液，此时气泡序号领先
		{
			n_liquid=1;
			n=n_trans[2*n_liquid_total];	//第1个液塞的左界面
			momentum=0;
			l_liquid=0;
			momentum=momentum+v_lright[n][2]*x_lright[n][2];
			l_liquid=l_liquid+x_lright[n][2];
			n=n+1;
			while (n<=length)	//液塞左界面右边控制体到坐标原点
			{
				momentum=momentum+v_lleft[n][2]*x_lleft[n][2];
				l_liquid=l_liquid+x_lleft[n][2];
				n=n+1;
			}
			n=1;
			while (n<=n_trans[1])	//坐标原点到右界面
			{
				momentum=momentum+v_lleft[n][2]*x_lleft[n][2];
				l_liquid=l_liquid+x_lleft[n][2];
				n=n+1;
			}
			n=n_trans[2*n_liquid_total];
			v_lright[n][2]=momentum/l_liquid;	//左界面新的速度
			n=n+1;
			while (n<=length)
			{
				v_lleft[n][2]=momentum/l_liquid;
				v_lright[n][2]=momentum/l_liquid;
				n=n+1; 
			}
			n=1;
			while (n<=n_trans[1]-1)
			{
				v_lleft[n][2]=momentum/l_liquid;
				v_lright[n][2]=momentum/l_liquid;
				n=n+1;
			}
			n=n_trans[1];
			v_lleft[n][2]=momentum/l_liquid;

			n_liquid=2;	//从第2个液塞开始，不会跨过原点
			while (n_liquid<=n_liquid_total)
			{
				n=n_trans[2*n_liquid-2];	//左界面
				momentum=0;
				l_liquid=0;
				momentum=momentum+v_lright[n][2]*x_lright[n][2];
				l_liquid=l_liquid+x_lright[n][2];
				n=n+1;
				while (n<=n_trans[2*n_liquid-1])	//左界面右边控制体到右界面
				{
					momentum=momentum+v_lleft[n][2]*x_lleft[n][2];
					l_liquid=l_liquid+x_lleft[n][2];
					n=n+1;
				}
				n=n_trans[2*n_liquid-2];
				v_lright[n][2]=momentum/l_liquid;	//左界面新的速度
				n=n+1;
				while (n<=n_trans[2*n_liquid-1]-1)
				{
					v_lleft[n][2]=momentum/l_liquid;
					v_lright[n][2]=momentum/l_liquid;
					n=n+1;
				}
				n=n_trans[2*n_liquid-1];
				v_lleft[n][2]=momentum/l_liquid;

				n_liquid=n_liquid+1;
			}
		}

		/*%%%%%%%%确定各气泡的状态--体积、质量、密度、温度、压力%%%%%%%%%*/
		//double c_vv[100];	//每个气泡的等容比热，J/kg-K
		//double c_pv[100];	//每个气泡的等压比热，J/kg-K
		double u_vb[100][3]; //每个气泡的内能，J/kg

		if (sort[1]==2||sort[1]==4)	//第1个控制体为气或气-液
		{
			n_bubble=1;
			while (n_bubble<=n_bubble_total)
			{
				//c_pv[n_bubble]=-1762920+399530*T_vb[n_bubble][1]-35873*pow(T_vb[n_bubble][1],2)+1605.24*pow(T_vb[n_bubble][1],3)-35.8187*pow(T_vb[n_bubble][1],4)+0.319537*pow(T_vb[n_bubble][1],5);
				//c_vv[n_bubble]=9711.08-615.358*T_vb[n_bubble][1]+41.7627*pow(T_vb[n_bubble][1],2)-1.2972*pow(T_vb[n_bubble][1],3)+0.0168863*pow(T_vb[n_bubble][1],4);

				if (n_bubble==1)	//第一个气泡单独考虑，因为它跨过坐标原点
				{
					V_bubble[n_bubble][2]=A_v*(x_v[n_trans[2*n_bubble-1]][2]+x_v[n_trans[2*n_bubble_total]][2]+(length-1+n_trans[2*n_bubble-1]-n_trans[2*n_bubble_total])*dx);
					m_bubble[n_bubble][2]=m_bubble[n_bubble][1]+Tr_v[n_bubble]*dt+Tr_vl[n_bubble]*dt;	//气泡质量=原质量+液膜蒸发+液体沸腾
					ro_vb[n_bubble][2]=m_bubble[n_bubble][2]/V_bubble[n_bubble][2];	//密度=质量/体积，而不是通过温度和压力查表得来
					/*气泡的能量方程*/
					//T_vb[n_bubble][2]=1/m_bubble[n_bubble][2]*(m_bubble[n_bubble][1]*T_vb[n_bubble][1]+dt/c_vv[n_bubble]*((Tr_v[n_bubble]+Tr_vl[n_bubble])*c_pv[n_bubble]*T_vb[n_bubble][1]-P_vb[n_bubble][1]/dt*(V_bubble[n_bubble][2]-V_bubble[n_bubble][1])));	//气泡内能变化=流进的焓-做功
					u_vb[n_bubble][1] = u_given_vt(V_bubble[n_bubble][1], m_bubble[n_bubble][1], T_vb[n_bubble][1]);//上一时刻的气泡内能
					u_vb[n_bubble][2] = u_vb[n_bubble][1] + 1.0 / m_bubble[n_bubble][1] * ((Tr_v[n_bubble] + Tr_vl[n_bubble])*dt*P_vb[n_bubble][1]* V_bubble[n_bubble][1]/ m_bubble[n_bubble][1]- P_vb[n_bubble][1] * (V_bubble[n_bubble][2] - V_bubble[n_bubble][1]));//气泡能量方程
					T_vb[n_bubble][2] = T_given_uv(u_vb[n_bubble][2], V_bubble[n_bubble][2], m_bubble[n_bubble][2], T_vb[n_bubble][1]);//新时刻的气泡温度
					P_vb[n_bubble][2]=mm*T_vb[n_bubble][2]/(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] -0.009136)-(6162* pow((1+0.1273*(1-sqrt(T_vb[n_bubble][2]/33.15))), 2))/(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] *(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] +0.009136));	//理想气体状态方程改为实际气体状态方程
					P_sat= p_sat_given_T(T_vb[n_bubble][2]);	//待确认
					if (P_vb[n_bubble][2]>P_sat)	//将过冷气体修正为饱和气体
					{
						P_vb[n_bubble][2]=P_sat;
						ro_vb[n_bubble][2]=ro_sat_given_T(T_vb[n_bubble][2]);	//将理想气体状态方程改为RKS的拟合式
						m_bubble[n_bubble][2]=ro_vb[n_bubble][2]*V_bubble[n_bubble][2];
					}
					T_sat_vb[n_bubble][2]= T_sat_given_p(P_vb[n_bubble][2]);//新增气泡饱和温度
					/*%对气泡内的各控制体进行初始化%*/
					n=1;	//从坐标原点到右界面
					while (n<=n_trans[2*n_bubble-1])
					{
						P_v[n][2]=P_vb[n_bubble][2];
						T_v[n][2]=T_vb[n_bubble][2];
						T_sat[n][2]=T_sat_vb[n_bubble][2];
						ro_v[n][2]=ro_vb[n_bubble][2];
						n=n+1;
					}
					n=n_trans[2*n_bubble_total];	//从左界面到坐标原点
					while (n<=length)
					{
						P_v[n][2]=P_vb[n_bubble][2];
						T_v[n][2]=T_vb[n_bubble][2];
						T_sat[n][2]=T_sat_vb[n_bubble][2];
						ro_v[n][2]=ro_vb[n_bubble][2];
						n=n+1;
					}
					
					Tr_v[n_bubble]=0;
					n=n_trans[2*n_bubble_total];	//从左界面到坐标原点
					while (n<=length)
					{
						P_sat= p_sat_given_T(T_w[n][1]);
						if (P_sat>=P_v[n][1])	//壁面温度对应的饱和压力高于气泡压力
						{
							Tr_v_sig[n]=1;	//蒸发
						}
						else
						{
							Tr_v_sig[n]=-1;	//冷凝
						}

						h_fg_v[n]= h_fg_given_T(T_v[n][1]);//为何这里有？
						
						Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//气泡蒸发、冷凝质量累加
						n=n+1;
					}
					n=1;
					while (n<=n_trans[2*n_bubble-1])	//从坐标原点到右界面
					{
						P_sat= p_sat_given_T(T_w[n][1]);
						if (P_sat>=P_v[n][1])	//壁面温度对应的饱和压力高于气泡压力
						{
							Tr_v_sig[n]=1;	//蒸发
						}
						else
						{
							Tr_v_sig[n]=-1;
						}

						h_fg_v[n]= h_fg_given_T(T_v[n][1]);//为何这里有？

						Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//气泡蒸发、冷凝质量累加
						n=n+1;
					}
				}
				else  //其他不跨过坐标原点的气泡
				{
					if (n_trans[2*n_bubble-1]==n_trans[2*n_bubble-2])	//气泡左右界面在同一个控制体
					{
						V_bubble[n_bubble][2]=A_v*x_v[n_trans[2*n_bubble-2]][2];
					}
					else  //左右界面不在同一控制体
					{
						V_bubble[n_bubble][2]=A_v*(x_v[n_trans[2*n_bubble-2]][2]+x_v[n_trans[2*n_bubble-1]][2]+(n_trans[2*n_bubble-1]-n_trans[2*n_bubble-2]-1)*dx);
					}
					m_bubble[n_bubble][2]=m_bubble[n_bubble][1]+Tr_v[n_bubble]*dt+Tr_vl[n_bubble]*dt;	//气泡新的质量
					ro_vb[n_bubble][2]=m_bubble[n_bubble][2]/V_bubble[n_bubble][2];	//气泡新的密度
					/*气泡的能量方程*/
					//T_vb[n_bubble][2]=1/m_bubble[n_bubble][2]*(m_bubble[n_bubble][1]*T_vb[n_bubble][1]+dt/c_vv[n_bubble]*((Tr_v[n_bubble]+Tr_vl[n_bubble])*c_pv[n_bubble]*T_vb[n_bubble][1]-P_vb[n_bubble][1]/dt*(V_bubble[n_bubble][2]-V_bubble[n_bubble][1])));	//气泡新的温度
					u_vb[n_bubble][1] = u_given_vt(V_bubble[n_bubble][1], m_bubble[n_bubble][1], T_vb[n_bubble][1]);//上一时刻的气泡内能
					u_vb[n_bubble][2] = u_vb[n_bubble][1] + 1.0 / m_bubble[n_bubble][1] * ((Tr_v[n_bubble] + Tr_vl[n_bubble])*dt*P_vb[n_bubble][1] * V_bubble[n_bubble][1] / m_bubble[n_bubble][1] - P_vb[n_bubble][1] * (V_bubble[n_bubble][2] - V_bubble[n_bubble][1]));//气泡能量方程
					T_vb[n_bubble][2] = T_given_uv(u_vb[n_bubble][2], V_bubble[n_bubble][2], m_bubble[n_bubble][2], T_vb[n_bubble][1]);//新时刻的气泡温度
					P_vb[n_bubble][2]=mm*T_vb[n_bubble][2]/(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] -0.009136)-(6162* pow((1+0.1273*(1-sqrt(T_vb[n_bubble][2]/33.15))), 2))/(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] *(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] +0.009136));	//理想气体状态方程改为实际气体状态方程
					P_sat= p_sat_given_T(T_vb[n_bubble][2]);	//气泡温度对应的饱和压力
					if (P_vb[n_bubble][2]>P_sat)	//过冷气泡，修正为饱和气泡
					{
						P_vb[n_bubble][2]=P_sat;
						ro_vb[n_bubble][2]= ro_sat_given_T(T_vb[n_bubble][2]);	//将理想气体状态方程改为RKS的拟合式
						m_bubble[n_bubble][2]=ro_vb[n_bubble][2]*V_bubble[n_bubble][2];
					}
					T_sat_vb[n_bubble][2]= T_sat_given_p(P_vb[n_bubble][2]);//新增气泡饱和温度
					/*%对气泡内的各控制体进行初始化%*/
					n=n_trans[2*n_bubble-2];	//左界面
					while (n<=n_trans[2*n_bubble-1])	//到右界面
					{
						P_v[n][2]=P_vb[n_bubble][2];
						T_v[n][2]=T_vb[n_bubble][2];
						T_sat[n][2] = T_sat_vb[n_bubble][2];
						ro_v[n][2]=ro_vb[n_bubble][2];
						n=n+1;
					}
				
					Tr_v[n_bubble]=0;
					n=n_trans[2*n_bubble-2];	//左界面
					while (n<=n_trans[2*n_bubble-1])	//到右界面
					{
						P_sat= p_sat_given_T(T_w[n][1]);
						if (P_sat>=P_v[n][1])	//壁面温度对应的饱和压力高于气泡压力
						{
							Tr_v_sig[n]=1;	//蒸发
						}
						else
						{
							Tr_v_sig[n]=-1;	//冷凝
						}
						h_fg_v[n]= h_fg_given_T(T_v[n][1]);//为何这里有？
						Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//气泡蒸发、冷凝质量累加
						n=n+1;
					}
				}
				
				/*%将气泡的状态值重新覆盖%*/
				P_vb[n_bubble][1]=P_vb[n_bubble][2];
				T_vb[n_bubble][1]=T_vb[n_bubble][2];
				T_sat_vb[n_bubble][1] = T_sat_vb[n_bubble][2];
				ro_vb[n_bubble][1]=ro_vb[n_bubble][2];
				V_bubble[n_bubble][1]=V_bubble[n_bubble][2];
				m_bubble[n_bubble][1]=m_bubble[n_bubble][2];
             	n_bubble=n_bubble+1;
			}
			/*液塞内压力线性分布*/
			n_bubble=1;
			while (n_bubble<=n_bubble_total)	//将气泡压力赋给液塞
			{
				P_l[n_trans[2*n_bubble-1]][2]=P_v[n_trans[2*n_bubble-1]][2];
				P_l[n_trans[2*n_bubble]][2]=P_v[n_trans[2*n_bubble]][2];
				n=n_trans[2*n_bubble-1]+1;	//液塞左界面右侧控制体
				while (n<=n_trans[2*n_bubble]-1)	//到液塞右界面左侧控制体
				{
					P_l[n][2]=P_l[n_trans[2*n_bubble]][2]-((n_trans[2*n_bubble]-n-0.5)*dx+x_lleft[n_trans[2*n_bubble]][2]/2)
						/((n_trans[2*n_bubble]-n_trans[2*n_bubble-1]-1)*dx+x_lleft[n_trans[2*n_bubble]][2]+x_lright[n_trans[2*n_bubble-1]][2])
						*(P_l[n_trans[2*n_bubble]][2]-P_l[n_trans[2*n_bubble-1]][2]);	//每个控制体内液体的压力表示在液体中心位置
					n=n+1;
				}
				n_bubble=n_bubble+1;
			}
		}
		else  //第1个控制体以液体开始，气泡序号领先液塞
		{
			n_bubble=1;
			while (n_bubble<=n_bubble_total)
			{
				//c_pv[n_bubble]=-1762920+399530*T_vb[n_bubble][1]-35873*pow(T_vb[n_bubble][1],2)+1605.24*pow(T_vb[n_bubble][1],3)-35.8187*pow(T_vb[n_bubble][1],4)+0.319537*pow(T_vb[n_bubble][1],5);
				//c_vv[n_bubble]=9711.08-615.358*T_vb[n_bubble][1]+41.7627*pow(T_vb[n_bubble][1],2)-1.2972*pow(T_vb[n_bubble][1],3)+0.0168863*pow(T_vb[n_bubble][1],4);
				if (n_trans[2*n_bubble]==n_trans[2*n_bubble-1])	//气泡新体积
				{
					V_bubble[n_bubble][2]=A_v*x_v[n_trans[2*n_bubble-1]][2];
				}
				else
				{
					V_bubble[n_bubble][2]=A_v*(x_v[n_trans[2*n_bubble-1]][2]+x_v[n_trans[2*n_bubble]][2]+(n_trans[2*n_bubble]-n_trans[2*n_bubble-1]-1)*dx);
				}
				m_bubble[n_bubble][2]=m_bubble[n_bubble][1]+Tr_v[n_bubble]*dt+Tr_vl[n_bubble]*dt;	//气泡新质量
				ro_vb[n_bubble][2]=m_bubble[n_bubble][2]/V_bubble[n_bubble][2];	//气泡新密度
				/*气泡的能量方程*/
				//T_vb[n_bubble][2]=1/m_bubble[n_bubble][2]*(m_bubble[n_bubble][1]*T_vb[n_bubble][1]+dt/c_vv[n_bubble]*((Tr_v[n_bubble]+Tr_vl[n_bubble])*c_pv[n_bubble]*T_vb[n_bubble][1]-P_vb[n_bubble][1]/dt*(V_bubble[n_bubble][2]-V_bubble[n_bubble][1])));	//原能量方程，废弃
				u_vb[n_bubble][1] = u_given_vt(V_bubble[n_bubble][1], m_bubble[n_bubble][1], T_vb[n_bubble][1]);//上一时刻的气泡内能
				u_vb[n_bubble][2] = u_vb[n_bubble][1] + 1.0 / m_bubble[n_bubble][1] * ((Tr_v[n_bubble] + Tr_vl[n_bubble])*dt*P_vb[n_bubble][1] * V_bubble[n_bubble][1] / m_bubble[n_bubble][1] - P_vb[n_bubble][1] * (V_bubble[n_bubble][2] - V_bubble[n_bubble][1]));//气泡能量方程
				T_vb[n_bubble][2] = T_given_uv(u_vb[n_bubble][2], V_bubble[n_bubble][2], m_bubble[n_bubble][2], T_vb[n_bubble][1]);//新时刻的气泡温度
				P_vb[n_bubble][2]=mm*T_vb[n_bubble][2]/(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] -0.009136)-(6162* pow((1+0.1273*(1-sqrt(T_vb[n_bubble][2]/33.15))), 2))/(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] *(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] +0.009136));	//理想气体状态方程改为实际气体状态方程
				P_sat= p_sat_given_T(T_vb[n_bubble][2]);	//气泡温度对应的饱和压力
				if (P_vb[n_bubble][2]>P_sat)	//过冷气泡被修正为饱和气泡
				{
					P_vb[n_bubble][2]=P_sat;
					ro_vb[n_bubble][2]= ro_sat_given_T(T_vb[n_bubble][2]);	//将理想气体状态方程改为RKS的拟合式
					m_bubble[n_bubble][2]=ro_vb[n_bubble][2]*V_bubble[n_bubble][2];
				}
				T_sat_vb[n_bubble][2]= T_sat_given_p(P_vb[n_bubble][2]);//新增气泡饱和温度
				/*%对气泡内的各控制体进行初始化%*/
				n=n_trans[2*n_bubble-1];	//气泡左界面
				while (n<=n_trans[2*n_bubble])	//到右界面
				{
					P_v[n][2]=P_vb[n_bubble][2];
					T_v[n][2]=T_vb[n_bubble][2];
					T_sat[n][2] = T_sat_vb[n_bubble][2];
					ro_v[n][2]=ro_vb[n_bubble][2];
					n=n+1;
				}
				
				Tr_v[n_bubble]=0;
				n=n_trans[2*n_bubble-1];
				while (n<=n_trans[2*n_bubble])
				{
					P_sat= p_sat_given_T(T_w[n][1]);	//壁面温度对应的饱和压力
					if (P_sat>=P_v[n][1])	//壁面温度对应的饱和压力大于气泡压力
					{
						Tr_v_sig[n]=1;	//蒸发
					}
					else
					{
						Tr_v_sig[n]=-1;	//冷凝
					}
					h_fg_v[n]= h_fg_given_T(T_v[n][1]);//为何这里有？
					Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//气泡内蒸发、冷凝质量
					n=n+1;
				}

				/*%将气泡的状态重新覆盖%*/
				P_vb[n_bubble][1]=P_vb[n_bubble][2];
				T_vb[n_bubble][1]=T_vb[n_bubble][2];
				T_sat_vb[n_bubble][1] = T_sat_vb[n_bubble][2];
				ro_vb[n_bubble][1]=ro_vb[n_bubble][2];
				V_bubble[n_bubble][1]=V_bubble[n_bubble][2];
				m_bubble[n_bubble][1]=m_bubble[n_bubble][2];

				n_bubble=n_bubble+1;
			}

			/*将气泡压力赋给液塞，液塞内部压力线性分布*/
			n_bubble=1;
			P_l[n_trans[2*n_bubble-1]][2]=P_v[n_trans[2*n_bubble-1]][2];	//液塞右界面压力
			P_l[n_trans[2*n_bubble_total]][2]=P_v[n_trans[2*n_bubble_total]][2];	//液塞左界面压力
			n=1;
			while (n<=n_trans[2*n_bubble-1]-1)	//坐标原点到右界面左侧控制体
			{
				P_l[n][2]=P_l[n_trans[2*n_bubble-1]][2]-((n_trans[2*n_bubble-1]-n-0.5)*dx+x_lleft[n_trans[2*n_bubble-1]][2]/2)
					/((n_trans[2*n_bubble-1]-n_trans[2*n_bubble_total]+length-1)*dx+x_lleft[n_trans[2*n_bubble-1]][2]/2+x_lright[n_trans[2*n_bubble_total]][2]/2)
					*(P_l[n_trans[2*n_bubble-1]][2]-P_l[n_trans[2*n_bubble_total]][2]);
				n=n+1;
			}
			n=n_trans[2*n_bubble_total]+1;	//液塞左界面
			while (n<=length)	//到坐标原点
			{
				P_l[n][2]=P_l[n_trans[2*n_bubble_total]][2]+((n-n_trans[2*n_bubble_total]-0.5)*dx+x_lright[n_trans[2*n_bubble_total]][2]/2)
					/((n_trans[2*n_bubble-1]-n_trans[2*n_bubble_total]+length-1)*dx+x_lleft[n_trans[2*n_bubble-1]][2]/2+x_lright[n_trans[2*n_bubble_total]][2]/2)
					*(P_l[n_trans[2*n_bubble-1]][2]-P_l[n_trans[2*n_bubble_total]][2]);
				n=n+1;
			}
			n_bubble=2;
			while (n_bubble<=n_bubble_total)
			{
				P_l[n_trans[2*n_bubble-1]][2]=P_v[n_trans[2*n_bubble-1]][2];	//液塞右界面压力
				P_l[n_trans[2*n_bubble-2]][2]=P_v[n_trans[2*n_bubble-2]][2];	//液塞左界面压力
				n=n_trans[2*n_bubble-2]+1;	//液塞左界面右边控制体
				while (n<=n_trans[2*n_bubble-1]-1)	//到液塞右界面左边控制体
				{
					P_l[n][2]=P_l[n_trans[2*n_bubble-1]][2]-((n_trans[2*n_bubble-1]-n-0.5)*dx+x_lleft[n_trans[2*n_bubble-1]][2]/2)
						/((n_trans[2*n_bubble-1]-n_trans[2*n_bubble-2]-1)*dx+x_lleft[n_trans[2*n_bubble-1]][2]/2+x_lright[n_trans[2*n_bubble-2]][2]/2)
						*(P_l[n_trans[2*n_bubble-1]][2]-P_l[n_trans[2*n_bubble-2]][2]);
					n=n+1;
				}
				n_bubble=n_bubble+1;
			}
		}

		/*****************气泡产生******************/
		double P_vb_new;	//新产生气泡的压力，Pa
		double T_vb_new;	//新产生气泡的温度，K
		double ro_vb_new;	//新产生气泡的密度，kg/m3
		double V_bubble_new;	//新产生气泡的体积，m3
		double m_bubble_new;	//新产生气泡的质量，kg
		double Tr_v_new;	//新产生气泡的蒸发、冷凝质量，kg
		int n_bubble_new;	//新产生气泡的序号

		/*%%第一个气泡产生点%%*/
		n=50;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)	//气泡生成点1的左右两边各有2个控制体的液体
		{
			if (generate_inter1>=generate_fre/dt)	//气泡生成点1的时间间隔达到要求
			{
				/*气泡生成点1过热度达到要求*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)//液体饱和温度待确定？？？
				{
					sort[n]=6;	
					x_lleft[n][2]=(dx-l_disappear)/2;	//气泡产生点在控制体正中间
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];	//新气泡温度=壁面温度
					T_sat[n][2] = T_v[n][2];	//新气泡饱和温度=气泡温度
					P_sat= p_sat_given_T(T_v[n][2]);	//新气泡饱和压力
					P_v[n][2]=P_sat;	//新气泡压力=饱和压力
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//将理想气体状态方程改为RKS的拟合式;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;	//质量根据密度和体积计算得来

					h_fg_v[n]= h_fg_given_T(T_v[n][1]);//为何这里有？

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);	//壁面温度对应的饱和压力
					if (P_sat>=P_v[n][2])	//壁面温度>=气泡温度，恒成立
					{
						Tr_v_sig[n]=1;	//蒸发
					}
					else  //不可能发生。。。。
					{
						Tr_v_sig[n]=-1;	//冷凝
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//气泡蒸发、冷凝质量累加
					
					generate_inter1=0;	//气泡产生之后，重置计时器
					bubble_gen_total=bubble_gen_total+1;	//新产生气泡数量+1

					/*重新判断新的气液交界面*/
					n=1;
					if (sort[n]==3||sort[n]==4)	//液-气或者气-液
					{
						transfer_point_number=1;
						n_trans[1]=1;
					}
					else if (sort[n]==8||sort[n]==6)	//液-气-液或者新气泡
					{
						transfer_point_number=2;
						n_trans[1]=1;
						n_trans[2]=1;
					}
					else  //液或者气
					{
						transfer_point_number=0;
					}

					n=2;
					while (n<=length)
					{
						if (sort[n]==3||sort[n]==4)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						if (sort[n]==8||sort[n]==6)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n;                  
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						n=n+1;
					}
					n_trans_total=transfer_point_number;                 /*%%记录交界面总数%%*/ 
					n_bubble_total=n_trans_total/2;                      /*%%气泡的数量%%*/
					n_liquid_total=n_trans_total/2;                      /*%%液塞的数量%%*/

					/*对气泡序号进行修正（由于气泡产生）*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//第1个控制体以液体开始
					{
						n_bubble=0;
						n_bubble_new=0;
						transfer_point_number = 1;	//从第1个气液界面开始，找到第1个新气泡
						while (n_bubble_new==0)	//直到找到新气泡为止
						{
							if ((transfer_point_number%2)!=0)	//奇数
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)	//该控制体有新气泡产生
								{
									n_bubble_new=n_bubble;	//第1个新气泡的序号，结束循环
									sort[n_trans[transfer_point_number]]=8;	//该控制体变成sort[n]=8
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					if (sort[1]==2||sort[1]==4)	//第1个控制体以气体开始
					{
						n_bubble=1;
						n_bubble_new=0;
						transfer_point_number = 1;
						while (n_bubble_new==0)	//直到找到新气泡
						{
							if ((transfer_point_number%2)==0)	//偶数
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)	//该控制体有气泡产生
								{
									n_bubble_new=n_bubble;	//第1个新气泡的序号，结束循环
									sort[n_trans[transfer_point_number]]=8;	//该控制体变成sort[n]=8
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					n_bubble=n_bubble_total;
					while (n_bubble>=n_bubble_new+1)	//第1个新气泡之后，每个气泡的序号增大1
					{
						P_vb[n_bubble][1]=P_vb[n_bubble-1][1];
						T_vb[n_bubble][1]=T_vb[n_bubble-1][1];
						T_sat_vb[n_bubble][1] = T_sat_vb[n_bubble - 1][1];
						ro_vb[n_bubble][1]=ro_vb[n_bubble-1][1];
						V_bubble[n_bubble][1]=V_bubble[n_bubble-1][1];
						m_bubble[n_bubble][1]=m_bubble[n_bubble-1][1];
						Tr_v[n_bubble]=Tr_v[n_bubble-1];
						n_bubble=n_bubble-1;
					}
					n_bubble=n_bubble_new;	//新气泡的参数
					P_vb[n_bubble][1]=P_vb_new;
					T_vb[n_bubble][1]=T_vb_new;
					T_sat_vb[n_bubble][1]= T_vb[n_bubble][1];
					ro_vb[n_bubble][1]=ro_vb_new;
					V_bubble[n_bubble][1]=V_bubble_new;
					m_bubble[n_bubble][1]=m_bubble_new;
					Tr_v[n_bubble]=Tr_v_new;
				}
			}
			else  //产生气泡的时间未到
			{
				generate_inter1=generate_inter1+1;	//时间步长数量+1
			}
		}
		else  //气泡生成点周围液体不够
		{
			generate_inter1=generate_inter1+1;	//步数计量器+1
		}

		/*%%第二个气泡产生点%%*/
		n=104;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)	//气泡产生点周围液体足够
		{
			if (generate_inter2>=generate_fre/dt)	//时间间隔达到要求
			{
				/*过热度达到要求*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)      /*%有新气泡产生%*/
				{
					sort[n]=6;	//有新气泡产生
					x_lleft[n][2]=(dx-l_disappear)/2;
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];
					T_sat[n][2]=T_v[n][2];
					P_sat= p_sat_given_T(T_v[n][2]);
					P_v[n][2]=P_sat;
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//将理想气体状态方程改为RKS的拟合式;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;

					h_fg_v[n]= h_fg_given_T(T_v[n][1]);//为何这里有？

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);
					if (P_sat>=P_v[n][2])	//这个等式恒成立
					{
						Tr_v_sig[n]=1;	//蒸发
					}
					else
					{
						Tr_v_sig[n]=-1;	//冷凝
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//新气泡蒸发、冷凝质量
					
					generate_inter2=0;	//重置计时器
					bubble_gen_total=bubble_gen_total+1;	//新产生的气泡总数量+1

					/*重新判断气液交界面位置*/
					n=1;
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=1;
						n_trans[1]=1;
					}
					else if (sort[n]==8||sort[n]==6)
					{
						transfer_point_number=2;
						n_trans[1]=1;
						n_trans[2]=1;
					}
					else
					{
						transfer_point_number=0;
					}

					n=2;
					while (n<=length)
					{
						if (sort[n]==3||sort[n]==4)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						if (sort[n]==8||sort[n]==6)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n;                  
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						n=n+1;
					}

					n_trans_total=transfer_point_number;	//气液界面总数
					n_bubble_total=n_trans_total/2;	//新的气泡总数
					n_liquid_total=n_trans_total/2;		//新的液塞数量

					/*对气泡序号进行修正（由于气泡产生）*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//第1个控制体以液体开始
					{
						transfer_point_number=1;
						n_bubble=0;
						n_bubble_new=0;
						while (n_bubble_new==0)	//直到找到新气泡
						{
							if ((transfer_point_number%2)!=0)
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)
								{
									n_bubble_new=n_bubble;
									sort[n_trans[transfer_point_number]]=8;
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					if (sort[1]==2||sort[1]==4)	//第1个控制体以气体开始
					{
						transfer_point_number=1;
						n_bubble=1;
						n_bubble_new=0;
						while (n_bubble_new==0)	//直到找到新气泡
						{
							if ((transfer_point_number%2)==0)
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)
								{
									n_bubble_new=n_bubble;
									sort[n_trans[transfer_point_number]]=8;
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					n_bubble=n_bubble_total;
					while (n_bubble>=n_bubble_new+1)	//新气泡之后的气泡序号+1
					{
						P_vb[n_bubble][1]=P_vb[n_bubble-1][1];
						T_vb[n_bubble][1]=T_vb[n_bubble-1][1];
						T_sat_vb[n_bubble][1] = T_sat_vb[n_bubble - 1][1];
						ro_vb[n_bubble][1]=ro_vb[n_bubble-1][1];
						V_bubble[n_bubble][1]=V_bubble[n_bubble-1][1];
						m_bubble[n_bubble][1]=m_bubble[n_bubble-1][1];
						Tr_v[n_bubble]=Tr_v[n_bubble-1];
						n_bubble=n_bubble-1;
					}
					n_bubble=n_bubble_new;
					P_vb[n_bubble][1]=P_vb_new;
					T_vb[n_bubble][1]=T_vb_new;
					T_sat_vb[n_bubble][1] = T_vb[n_bubble][1];
					ro_vb[n_bubble][1]=ro_vb_new;
					V_bubble[n_bubble][1]=V_bubble_new;
					m_bubble[n_bubble][1]=m_bubble_new;
					Tr_v[n_bubble]=Tr_v_new;
				}
			}
			else  //时间间隔不够
			{
				generate_inter2=generate_inter2+1;	//时间步长增加1
			}
		}
		else  //气泡产生点周围液体不够
		{
			generate_inter2=generate_inter2+1;	//步数计量器+1
		}

		/*%%第三个气泡产生点%%*/
		n=206;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)	//气泡产生点左右液体足够
		{
			if (generate_inter3>=generate_fre/dt)	//时间间隔达到要求
			{
				/*过热度达到要求*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)      /*%有新气泡产生%*/
				{
					sort[n]=6;
					x_lleft[n][2]=(dx-l_disappear)/2;
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];	//新气泡温度=壁面温度
					T_sat[n][2]=T_v[n][2];
					P_sat= p_sat_given_T(T_v[n][2]);
					P_v[n][2]=P_sat;
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//将理想气体状态方程改为RKS的拟合式;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;

					h_fg_v[n] = h_fg_given_T(T_v[n][1]);//为何这里有？

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);
					if (P_sat>=P_v[n][2])	//等式恒成立
					{
						Tr_v_sig[n]=1;	//蒸发
					}
					else
					{
						Tr_v_sig[n]=-1;
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//新气泡蒸发、冷凝质量
					
					generate_inter3=0;	//重置计时器
					bubble_gen_total=bubble_gen_total+1;	//新产生的气泡数量+1

					/*重新判断气液交界面位置*/
					n=1;
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=1;
						n_trans[1]=1;
					}
					else if (sort[n]==8||sort[n]==6)
					{
						transfer_point_number=2;
						n_trans[1]=1;
						n_trans[2]=1;
					}
					else
					{
						transfer_point_number=0;
					}

					n=2;
					while (n<=length)
					{
						if (sort[n]==3||sort[n]==4)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						if (sort[n]==8||sort[n]==6)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n;                  
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						n=n+1;
					}

					n_trans_total=transfer_point_number;                 /*%%记录交界面总数%%*/ 
					n_bubble_total=n_trans_total/2;                      /*%%气泡的数量%%*/
					n_liquid_total=n_trans_total/2;                      /*%%液塞的数量%%*/

					/*对气泡序号进行修正（由于气泡产生）*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//第1个控制体以液开始
					{
						transfer_point_number=1;
						n_bubble=0;
						n_bubble_new=0;
						while (n_bubble_new==0)	//直到找到新气泡
						{
							if ((transfer_point_number%2)!=0)
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)
								{
									n_bubble_new=n_bubble;
									sort[n_trans[transfer_point_number]]=8;
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					if (sort[1]==2||sort[1]==4)	//第1个控制体以气开始
					{
						transfer_point_number=1;
						n_bubble=1;
						n_bubble_new=0;
						while (n_bubble_new==0)	//直到找到新气泡序号
						{
							if ((transfer_point_number%2)==0)
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)
								{
									n_bubble_new=n_bubble;
									sort[n_trans[transfer_point_number]]=8;
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					n_bubble=n_bubble_total;
					while (n_bubble>=n_bubble_new+1)	//新气泡之后的气泡序号+1
					{
						P_vb[n_bubble][1]=P_vb[n_bubble-1][1];
						T_vb[n_bubble][1]=T_vb[n_bubble-1][1];
						T_sat_vb[n_bubble][1] = T_sat_vb[n_bubble - 1][1];
						ro_vb[n_bubble][1]=ro_vb[n_bubble-1][1];
						V_bubble[n_bubble][1]=V_bubble[n_bubble-1][1];
						m_bubble[n_bubble][1]=m_bubble[n_bubble-1][1];
						Tr_v[n_bubble]=Tr_v[n_bubble-1];
						n_bubble=n_bubble-1;
					}
					n_bubble=n_bubble_new;
					P_vb[n_bubble][1]=P_vb_new;
					T_vb[n_bubble][1]=T_vb_new;
					T_sat_vb[n_bubble][1] = T_vb[n_bubble][1];
					ro_vb[n_bubble][1]=ro_vb_new;
					V_bubble[n_bubble][1]=V_bubble_new;
					m_bubble[n_bubble][1]=m_bubble_new;
					Tr_v[n_bubble]=Tr_v_new;
				}
			}
			else  //时间间隔不够
			{
				generate_inter3=generate_inter3+1;	//时间步长数量+1
			}
		}
		else  //气泡产生点左右液体不够
		{
			generate_inter3=generate_inter3+1;	//步数计量器+1
		}

		/*%%第四个气泡产生点%%*/
		n=260;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)	//气泡产生点左右液体足够
		{
			if (generate_inter4>=generate_fre/dt)	//时间间隔满足要求
			{
				/*过热度满足要求*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)      /*%有新气泡产生%*/
				{
					sort[n]=6;
					x_lleft[n][2]=(dx-l_disappear)/2;
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];	//新气泡温度=壁面温度
					T_sat[n][2] = T_v[n][2];
					P_sat= p_sat_given_T(T_v[n][2]);
					P_v[n][2]=P_sat;
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//将理想气体状态方程改为RKS的拟合式;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;

					h_fg_v[n] = h_fg_given_T(T_v[n][1]);;//为何这里有？

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);
					if (P_sat>=P_v[n][2])	//等式恒成立
					{
						Tr_v_sig[n]=1;	//蒸发
					}
					else
					{
						Tr_v_sig[n]=-1;
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//气泡蒸发、冷凝质量
					
					generate_inter4=0;	//重置步数计量器
					bubble_gen_total=bubble_gen_total+1;	//新产生气泡数量+1

					/*重新判断气液交界面位置*/
					n=1;
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=1;
						n_trans[1]=1;
					}
					else if (sort[n]==8||sort[n]==6)
					{
						transfer_point_number=2;
						n_trans[1]=1;
						n_trans[2]=1;
					}
					else
					{
						transfer_point_number=0;
					}

					n=2;
					while (n<=length)
					{
						if (sort[n]==3||sort[n]==4)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						if (sort[n]==8||sort[n]==6)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n;                  
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						n=n+1;
					}

					n_trans_total=transfer_point_number;                 /*%%记录交界面总数%%*/ 
					n_bubble_total=n_trans_total/2;                      /*%%气泡的数量%%*/
					n_liquid_total=n_trans_total/2;                      /*%%液塞的数量%%*/

					/*对气泡序号进行修正（由于气泡产生）*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//第1个控制体以液开始
					{
						transfer_point_number=1;
						n_bubble=0;
						n_bubble_new=0;
						while (n_bubble_new==0)	//直至找到新气泡的序号
						{
							if ((transfer_point_number%2)!=0)	//奇数
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)
								{
									n_bubble_new=n_bubble;
									sort[n_trans[transfer_point_number]]=8;
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					if (sort[1]==2||sort[1]==4)	//第1个控制体以气开始
					{
						transfer_point_number=1;
						n_bubble=1;
						n_bubble_new=0;
						while (n_bubble_new==0)	//直至找到新气泡位置
						{
							if ((transfer_point_number%2)==0)	//偶数
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)
								{
									n_bubble_new=n_bubble;
									sort[n_trans[transfer_point_number]]=8;
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					n_bubble=n_bubble_total;
					while (n_bubble>=n_bubble_new+1)	//新气泡之后的气泡序号+1
					{
						P_vb[n_bubble][1]=P_vb[n_bubble-1][1];
						T_vb[n_bubble][1]=T_vb[n_bubble-1][1];
						T_sat_vb[n_bubble][1] = T_sat_vb[n_bubble - 1][1];
						ro_vb[n_bubble][1]=ro_vb[n_bubble-1][1];
						V_bubble[n_bubble][1]=V_bubble[n_bubble-1][1];
						m_bubble[n_bubble][1]=m_bubble[n_bubble-1][1];
						Tr_v[n_bubble]=Tr_v[n_bubble-1];
						n_bubble=n_bubble-1;
					}
					n_bubble=n_bubble_new;
					P_vb[n_bubble][1]=P_vb_new;
					T_vb[n_bubble][1]=T_vb_new;
					T_sat_vb[n_bubble][1] = T_vb[n_bubble][1];
					ro_vb[n_bubble][1]=ro_vb_new;
					V_bubble[n_bubble][1]=V_bubble_new;
					m_bubble[n_bubble][1]=m_bubble_new;
					Tr_v[n_bubble]=Tr_v_new;
				}
			}
			else  //时间间隔不满足要求
			{
				generate_inter4=generate_inter4+1;
			}
		}
		else  //气泡产生点左右液体不够
		{
			generate_inter4=generate_inter4+1;	//步数计量器+1
		}

		/*%%第五个气泡产生点%%*/
		n=470;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)
		{
			if (generate_inter5>=generate_fre/dt)
			{
				/*过热度满足要求*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)      /*%有新气泡产生%*/
				{
					sort[n]=6;
					x_lleft[n][2]=(dx-l_disappear)/2;
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];	//新气泡温度=壁面温度
					T_sat[n][2] = T_v[n][2];
					P_sat= p_sat_given_T(T_v[n][2]);
					P_v[n][2]=P_sat;
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//将理想气体状态方程改为RKS的拟合式;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;

					h_fg_v[n] = h_fg_given_T(T_v[n][1]);;//为何这里有？

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);
					if (P_sat>=P_v[n][2])	//等式恒成立
					{
						Tr_v_sig[n]=1;	//蒸发
					}
					else
					{
						Tr_v_sig[n]=-1;
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//新气泡蒸发、冷凝质量
					
					generate_inter5=0;	//重置步数计量器
					bubble_gen_total=bubble_gen_total+1;	//新产生气泡数量+1

					/*重新判断气液交界面位置*/
					n=1;
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=1;
						n_trans[1]=1;
					}
					else if (sort[n]==8||sort[n]==6)
					{
						transfer_point_number=2;
						n_trans[1]=1;
						n_trans[2]=1;
					}
					else
					{
						transfer_point_number=0;
					}

					n=2;
					while (n<=length)
					{
						if (sort[n]==3||sort[n]==4)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						if (sort[n]==8||sort[n]==6)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n;                  
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						n=n+1;
					}

					n_trans_total=transfer_point_number;                 /*%%记录交界面总数%%*/ 
					n_bubble_total=n_trans_total/2;                      /*%%气泡的数量%%*/
					n_liquid_total=n_trans_total/2;                      /*%%液塞的数量%%*/

					/*对气泡序号进行修正（由于气泡产生）*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//第1个控制体以液开始
					{
						transfer_point_number=1;
						n_bubble=0;
						n_bubble_new=0;
						while (n_bubble_new==0)	//直至找到新气泡序号
						{
							if ((transfer_point_number%2)!=0)	//奇数
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)
								{
									n_bubble_new=n_bubble;
									sort[n_trans[transfer_point_number]]=8;
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					if (sort[1]==2||sort[1]==4)	//第1个控制体以气开始
					{
						transfer_point_number=1;
						n_bubble=1;
						n_bubble_new=0;
						while (n_bubble_new==0)	//直至找到新气泡序号
						{
							if ((transfer_point_number%2)==0)	//偶数
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)
								{
									n_bubble_new=n_bubble;
									sort[n_trans[transfer_point_number]]=8;
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					n_bubble=n_bubble_total;
					while (n_bubble>=n_bubble_new+1)	//新气泡之后的气泡序号+1
					{
						P_vb[n_bubble][1]=P_vb[n_bubble-1][1];
						T_vb[n_bubble][1]=T_vb[n_bubble-1][1];
						T_sat_vb[n_bubble][1] = T_sat_vb[n_bubble - 1][1];
						ro_vb[n_bubble][1]=ro_vb[n_bubble-1][1];
						V_bubble[n_bubble][1]=V_bubble[n_bubble-1][1];
						m_bubble[n_bubble][1]=m_bubble[n_bubble-1][1];
						Tr_v[n_bubble]=Tr_v[n_bubble-1];
						n_bubble=n_bubble-1;
					}
					n_bubble=n_bubble_new;
					P_vb[n_bubble][1]=P_vb_new;
					T_vb[n_bubble][1]=T_vb_new;
					T_sat_vb[n_bubble][1] = T_vb[n_bubble][1];
					ro_vb[n_bubble][1]=ro_vb_new;
					V_bubble[n_bubble][1]=V_bubble_new;
					m_bubble[n_bubble][1]=m_bubble_new;
					Tr_v[n_bubble]=Tr_v_new;
				}

			}
			else  //时间间隔不够
			{
				generate_inter5=generate_inter5+1;	//步数计量器+1
			}
		}
		else  //气泡产生点左右液体不够
		{
			generate_inter5=generate_inter5+1;	//步数计量器+1
		}

		/*%%第六个气泡产生点%%*/
		n=524;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)
		{
			if (generate_inter6>=generate_fre/dt)
			{
				/*过热度满足要求*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)      /*%有新气泡产生%*/
				{
					sort[n]=6;
					x_lleft[n][2]=(dx-l_disappear)/2;
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];	//新气泡温度=壁面温度
					T_sat[n][2] = T_v[n][2];
					P_sat= p_sat_given_T(T_v[n][2]);
					P_v[n][2]=P_sat;
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//将理想气体状态方程改为RKS的拟合式;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;

					h_fg_v[n] = h_fg_given_T(T_v[n][1]);//为何这里有？

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);
					if (P_sat>=P_v[n][2])	//等式恒成立
					{
						Tr_v_sig[n]=1;	//蒸发
					}
					else
					{
						Tr_v_sig[n]=-1;
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//新气泡蒸发、冷凝质量

					generate_inter6=0;	//重置步数计量器
					bubble_gen_total=bubble_gen_total+1;	//新气泡数量+1

					/*重新判断气液交界面位置*/
					n=1;
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=1;
						n_trans[1]=1;
					}
					else if (sort[n]==8||sort[n]==6)
					{
						transfer_point_number=2;
						n_trans[1]=1;
						n_trans[2]=1;
					}
					else
					{
						transfer_point_number=0;
					}

					n=2;
					while (n<=length)
					{
						if (sort[n]==3||sort[n]==4)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						if (sort[n]==8||sort[n]==6)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n;                  
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						n=n+1;
					}

					n_trans_total=transfer_point_number;                 /*%%记录交界面总数%%*/ 
					n_bubble_total=n_trans_total/2;                      /*%%气泡的数量%%*/
					n_liquid_total=n_trans_total/2;                      /*%%液塞的数量%%*/

					/*对气泡序号进行修正（由于气泡产生）*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//第1个控制体以液开始
					{
						transfer_point_number=1;
						n_bubble=0;
						n_bubble_new=0;
						while (n_bubble_new==0)	//直至找到新气泡序号
						{
							if ((transfer_point_number%2)!=0)	//奇数
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)
								{
									n_bubble_new=n_bubble;
									sort[n_trans[transfer_point_number]]=8;
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					if (sort[1]==2||sort[1]==4)	//第1个控制体以气开始
					{
						transfer_point_number=1;
						n_bubble=1;
						n_bubble_new=0;
						while (n_bubble_new==0)	//直至找到新气泡序号
						{
							if ((transfer_point_number%2)==0)	//偶数
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)
								{
									n_bubble_new=n_bubble;
									sort[n_trans[transfer_point_number]]=8;
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					n_bubble=n_bubble_total;
					while (n_bubble>=n_bubble_new+1)	//新气泡之后的气泡序号+1
					{
						P_vb[n_bubble][1]=P_vb[n_bubble-1][1];
						T_vb[n_bubble][1]=T_vb[n_bubble-1][1];
						T_sat_vb[n_bubble][1] = T_sat_vb[n_bubble - 1][1];
						ro_vb[n_bubble][1]=ro_vb[n_bubble-1][1];
						V_bubble[n_bubble][1]=V_bubble[n_bubble-1][1];
						m_bubble[n_bubble][1]=m_bubble[n_bubble-1][1];
						Tr_v[n_bubble]=Tr_v[n_bubble-1];
						n_bubble=n_bubble-1;
					}
					n_bubble=n_bubble_new;
					P_vb[n_bubble][1]=P_vb_new;
					T_vb[n_bubble][1]=T_vb_new;
					T_sat_vb[n_bubble][1] = T_vb[n_bubble][1];
					ro_vb[n_bubble][1]=ro_vb_new;
					V_bubble[n_bubble][1]=V_bubble_new;
					m_bubble[n_bubble][1]=m_bubble_new;
					Tr_v[n_bubble]=Tr_v_new;
				}

			}
			else  //时间间隔不够
			{
				generate_inter6=generate_inter6+1;
			}
		}
		else  //气泡产生点左右液体不够
		{
			generate_inter6=generate_inter6+1;
		}

		/*%%第七个气泡产生点%%*/
		n=626;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)
		{
			if (generate_inter7>=generate_fre/dt)
			{
				/*过热度满足要求*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)      /*%有新气泡产生%*/
				{
					sort[n]=6;
					x_lleft[n][2]=(dx-l_disappear)/2;
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];	//新气泡温度=壁面温度
					T_sat[n][2] = T_v[n][2];
					P_sat= p_sat_given_T(T_v[n][2]);
					P_v[n][2]=P_sat;
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//将理想气体状态方程改为RKS的拟合式;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;

					h_fg_v[n] = h_fg_given_T(T_v[n][1]);//为何这里有？

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);
					if (P_sat>=P_v[n][2])	//等式恒成立
					{
						Tr_v_sig[n]=1;	//蒸发
					}
					else
					{
						Tr_v_sig[n]=-1;
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//新气泡蒸发、冷凝质量
					
					generate_inter7=0;	//重置步数计量器
					bubble_gen_total=bubble_gen_total+1;	//新产生气泡数量+1

					/*重新判断气液交界面位置*/
					n=1;
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=1;
						n_trans[1]=1;
					}
					else if (sort[n]==8||sort[n]==6)
					{
						transfer_point_number=2;
						n_trans[1]=1;
						n_trans[2]=1;
					}
					else
					{
						transfer_point_number=0;
					}

					n=2;
					while (n<=length)
					{
						if (sort[n]==3||sort[n]==4)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						if (sort[n]==8||sort[n]==6)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n;                  
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						n=n+1;
					}

					n_trans_total=transfer_point_number;                 /*%%记录交界面总数%%*/ 
					n_bubble_total=n_trans_total/2;                      /*%%气泡的数量%%*/
					n_liquid_total=n_trans_total/2;                      /*%%液塞的数量%%*/

					/*对气泡序号重新修正（由于气泡产生）*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//第1个控制体以液开始
					{
						transfer_point_number=1;
						n_bubble=0;
						n_bubble_new=0;
						while (n_bubble_new==0)	//直至找到新气泡序号
						{
							if ((transfer_point_number%2)!=0)	//奇数
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)
								{
									n_bubble_new=n_bubble;
									sort[n_trans[transfer_point_number]]=8;
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					if (sort[1]==2||sort[1]==4)	//第1个控制体以气开始
					{
						transfer_point_number=1;
						n_bubble=1;
						n_bubble_new=0;
						while (n_bubble_new==0)	//直至找到新气泡序号
						{
							if ((transfer_point_number%2)==0)	//偶数
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)
								{
									n_bubble_new=n_bubble;
									sort[n_trans[transfer_point_number]]=8;
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					n_bubble=n_bubble_total;
					while (n_bubble>=n_bubble_new+1)	//新气泡之后的气泡序号+1
					{
						P_vb[n_bubble][1]=P_vb[n_bubble-1][1];
						T_vb[n_bubble][1]=T_vb[n_bubble-1][1];
						T_sat_vb[n_bubble][1] = T_sat_vb[n_bubble - 1][1];
						ro_vb[n_bubble][1]=ro_vb[n_bubble-1][1];
						V_bubble[n_bubble][1]=V_bubble[n_bubble-1][1];
						m_bubble[n_bubble][1]=m_bubble[n_bubble-1][1];
						Tr_v[n_bubble]=Tr_v[n_bubble-1];
						n_bubble=n_bubble-1;
					}
					n_bubble=n_bubble_new;
					P_vb[n_bubble][1]=P_vb_new;
					T_vb[n_bubble][1]=T_vb_new;
					T_sat_vb[n_bubble][1] = T_vb[n_bubble][1];
					ro_vb[n_bubble][1]=ro_vb_new;
					V_bubble[n_bubble][1]=V_bubble_new;
					m_bubble[n_bubble][1]=m_bubble_new;
					Tr_v[n_bubble]=Tr_v_new;
				}

			}
			else  //时间间隔不够
			{
				generate_inter7=generate_inter7+1;
			}
		}
		else  //气泡产生点左右液体不够
		{
			generate_inter7=generate_inter7+1;
		}

		/*%%第八个气泡产生点%%*/
		n=680;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)
		{
			if (generate_inter8>=generate_fre/dt)
			{
				/*过热度满足条件*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)      /*%有新气泡产生%*/
				{
					sort[n]=6;
					x_lleft[n][2]=(dx-l_disappear)/2;
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];	//新气泡温度=壁面温度
					T_sat[n][2] = T_v[n][2];
					P_sat= p_sat_given_T(T_v[n][2]);
					P_v[n][2]=P_sat;
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//将理想气体状态方程改为RKS的拟合式;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;

					h_fg_v[n] = h_fg_given_T(T_v[n][1]);//为何这里有？

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);
					if (P_sat>=P_v[n][2])	//等式恒成立
					{
						Tr_v_sig[n]=1;	//蒸发
					}
					else
					{
						Tr_v_sig[n]=-1;
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//新气泡蒸发、冷凝质量
					
					generate_inter8=0;	//重置步数计量器
					bubble_gen_total=bubble_gen_total+1;	//新产生气泡数量+1

					/*重新判断气液交界面位置*/
					n=1;
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=1;
						n_trans[1]=1;
					}
					else if (sort[n]==8||sort[n]==6)
					{
						transfer_point_number=2;
						n_trans[1]=1;
						n_trans[2]=1;
					}
					else
					{
						transfer_point_number=0;
					}

					n=2;
					while (n<=length)
					{
						if (sort[n]==3||sort[n]==4)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						if (sort[n]==8||sort[n]==6)
						{
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n;                  
							transfer_point_number=transfer_point_number+1;
							n_trans[transfer_point_number]=n; 
						}
						n=n+1;
					}

					n_trans_total=transfer_point_number;                 /*%%记录交界面总数%%*/ 
					n_bubble_total=n_trans_total/2;                      /*%%气泡的数量%%*/
					n_liquid_total=n_trans_total/2;                      /*%%液塞的数量%%*/

					/*对气泡序号进行修正（由于气泡产生）*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//第1个控制体以液开始
					{
						transfer_point_number=1;
						n_bubble=0;
						n_bubble_new=0;
						while (n_bubble_new==0)	//直至找到新气泡序号
						{
							if ((transfer_point_number%2)!=0)	//奇数
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)
								{
									n_bubble_new=n_bubble;
									sort[n_trans[transfer_point_number]]=8;
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					if (sort[1]==2||sort[1]==4)	//第1个控制体以气开始
					{
						transfer_point_number=1;
						n_bubble=1;
						n_bubble_new=0;
						while (n_bubble_new==0)	//直至找到第1个气泡
						{
							if ((transfer_point_number%2)==0)	//偶数
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)
								{
									n_bubble_new=n_bubble;
									sort[n_trans[transfer_point_number]]=8;
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					n_bubble=n_bubble_total;
					while (n_bubble>=n_bubble_new+1)	//新气泡之后的气泡序号+1
					{
						P_vb[n_bubble][1]=P_vb[n_bubble-1][1];
						T_vb[n_bubble][1]=T_vb[n_bubble-1][1];
						T_sat_vb[n_bubble][1] = T_sat_vb[n_bubble - 1][1];
						ro_vb[n_bubble][1]=ro_vb[n_bubble-1][1];
						V_bubble[n_bubble][1]=V_bubble[n_bubble-1][1];
						m_bubble[n_bubble][1]=m_bubble[n_bubble-1][1];
						Tr_v[n_bubble]=Tr_v[n_bubble-1];
						n_bubble=n_bubble-1;
					}
					n_bubble=n_bubble_new;
					P_vb[n_bubble][1]=P_vb_new;
					T_vb[n_bubble][1]=T_vb_new;
					T_sat_vb[n_bubble][1] = T_vb[n_bubble][1];
					ro_vb[n_bubble][1]=ro_vb_new;
					V_bubble[n_bubble][1]=V_bubble_new;
					m_bubble[n_bubble][1]=m_bubble_new;
					Tr_v[n_bubble]=Tr_v_new;
				}

			}
			else  //时间间隔不够
			{
				generate_inter8=generate_inter8+1;
			}
		}
		else  //气泡产生点左右液体不够
		{
			generate_inter8=generate_inter8+1;
		}

		/**********************对旧值进行覆盖********************/
		n=1;
		/*新时刻控制体内气液长度，速度，温度，压力，密度*/
		while (n<=length)
		{
			x_lleft[n][1]=x_lleft[n][2];
			x_lright[n][1]=x_lright[n][2];
			x_v[n][1]=x_v[n][2];
			v_lleft[n][1]=v_lleft[n][2];
			v_lright[n][1]=v_lright[n][2];
			T_lleft[n][1]=T_lleft[n][2];
			T_lright[n][1]=T_lright[n][2];
			T_w[n][1]=T_w[n][2];
			P_v[n][1]=P_v[n][2];
			ro_v[n][1]=ro_v[n][2];
			T_v[n][1]=T_v[n][2];
			T_sat[n][1] = T_sat[n][2];
			n=n+1;
		}
		/*新时刻和温度有关的物性*/
		for (n=1;n<=length;n=n+1)
		{
			if (sort[n]==1)	//液
			{
				mu_left[n]=mu_l_given_T(T_lleft[n][1]);
				mu_right[n]=mu_l_given_T(T_lright[n][1]);
				c_pl_left[n]=cp_l_given_T(T_lleft[n][1]);
				c_pl_right[n]=cp_l_given_T(T_lright[n][1]);
				lamt_left[n]=lamt_l_given_T(T_lleft[n][1]);
				lamt_right[n]=lamt_l_given_T(T_lright[n][1]);
				h_fg_lleft[n]= h_fg_given_T(T_lleft[n][1]);
				h_fg_lright[n]= h_fg_given_T(T_lright[n][1]);
				h_fg_v[n]=0;
			}

			if (sort[n]==2)	//气
			{
				mu_left[n]=0;
				mu_right[n]=0;
				c_pl_left[n]=0;
				c_pl_right[n]=0;
				lamt_left[n]=0;
				lamt_right[n]=0;
				h_fg_lleft[n]=0;
				h_fg_lright[n]=0;
				h_fg_v[n] = h_fg_given_T(T_v[n][1]);
			}

			if (sort[n]==3)	//液-气
			{
				mu_left[n]=mu_l_given_T(T_lleft[n][1]);
				mu_right[n]=0;
				c_pl_left[n]=cp_l_given_T(T_lleft[n][1]);
				c_pl_right[n]=0;
				lamt_left[n]=lamt_l_given_T(T_lleft[n][1]);
				lamt_right[n]=0;
				h_fg_lleft[n]= h_fg_given_T(T_lleft[n][1]);
				h_fg_lright[n]=0;
				h_fg_v[n]= h_fg_given_T(T_v[n][1]);
			}

			if (sort[n]==4)	//气-液
			{
				mu_left[n]=0;
				mu_right[n]=mu_l_given_T(T_lright[n][1]);
				c_pl_left[n]=0;
				c_pl_right[n]=cp_l_given_T(T_lright[n][1]);
				lamt_left[n]=0;
				lamt_right[n]=lamt_l_given_T(T_lright[n][1]);
				h_fg_lleft[n]=0;
				h_fg_lright[n]= h_fg_given_T(T_lright[n][1]);
				h_fg_v[n]= h_fg_given_T(T_v[n][1]);
			}
			if (sort[n]==8)	//液-气-液
			{
				mu_left[n]=mu_l_given_T(T_lleft[n][1]);
				mu_right[n]=mu_l_given_T(T_lright[n][1]);
				c_pl_left[n]=cp_l_given_T(T_lleft[n][1]);
				c_pl_right[n]=cp_l_given_T(T_lright[n][1]);
				lamt_left[n]=lamt_l_given_T(T_lleft[n][1]);
				lamt_right[n]=lamt_l_given_T(T_lright[n][1]);
				h_fg_lleft[n]= h_fg_given_T(T_lleft[n][1]);
				h_fg_lright[n]= h_fg_given_T(T_lright[n][1]);
				h_fg_v[n] = h_fg_given_T(T_v[n][1]);
			}
		}
		/*新时刻的液体单相换热系数*/
		for(n=1;n<=length;n=n+1)
		{
			if (sort[n]==1||sort[n]==3||sort[n]==8)	//以液体开始，左侧液体的换热系数
			{
				Re_lleft[n]=d*ro_l*fabs(v_lleft[n][1])/mu_left[n];
				Pr = Pr_l_given_T(T_lleft[n][1]);
				if (Re_lleft[n]<=2300)
				{
					h_lleft[n]=4.364*lamt_left[n]/d;
				}
				else if (Re_lleft[n]<=10000)
				{
					f=pow(1.82*log10(Re_lleft[n])-1.64,-2);
					h_lleft[n]=lamt_left[n]/d*(f/8)*(Re_lleft[n]-1000)*Pr/(1+12.7*pow(f/8,0.5)*(pow(Pr,2/3)-1));
				}
				else
				{
					if (T_w[n][1]>T_lleft[n][1])
					{
						h_lleft[n]=0.023*lamt_left[n]/d*pow(Re_lleft[n],0.8)*pow(Pr,0.4);
					}
					else
					{
						h_lleft[n]=0.023*lamt_left[n]/d*pow(Re_lleft[n],0.8)*pow(Pr,0.3);
					}
				}
			}
			else  //以气体开始
			{
				Re_lleft[n]=0;
				h_lleft[n]=0;
			}
		}
		for(n=1;n<=length;n=n+1)
		{
			if (sort[n]==1||sort[n]==4||sort[n]==8)	//右侧液体的换热系数
			{
				Re_lright[n]=d*ro_l*fabs(v_lright[n][1])/mu_right[n];
				Pr = Pr_l_given_T(T_lright[n][1]);
				if (Re_lright[n]<=2300)
				{
					h_lright[n]=4.364*lamt_right[n]/d;
				}
				else if (Re_lright[n]<=10000)
				{
					f=pow(1.82*log10(Re_lright[n])-1.64,-2);
					h_lright[n]=lamt_right[n]/d*(f/8)*(Re_lright[n]-1000)*Pr/(1+12.7*pow(f/8,0.5)*(pow(Pr,2/3)-1));
				}
				else
				{
					if (T_w[n][1]>T_lright[n][1])
					{
						h_lright[n]=0.023*lamt_right[n]/d*pow(Re_lright[n],0.8)*pow(Pr,0.4);
					}
					else
					{
						h_lright[n]=0.023*lamt_right[n]/d*pow(Re_lright[n],0.8)*pow(Pr,0.3);
					}
				}
			}
			else  //右侧没有液体
			{
				Re_lright[n]=0;
				h_lright[n]=0;
			}
		}
		/*新时刻的液体沸腾传热系数*/
		for (n=1;n<=99;n=n+1)	//每个液塞（最多99个）初始化
		{
			Tr_vl_lleft[n]=0;
			Tr_vl_lright[n]=0;
			Tr_vl[n]=0;
			dl_vl_left[n]=0;
			dl_vl_right[n]=0;
		}
		for (n=1;n<=length;n=n+1)
		{
			Tr_vl_right_each[n]=0;
			Tr_vl_left_each[n]=0;
			h_b[n]=0; //沸腾换热系数初始化
		}
		/****根据液塞左右界面的位置来计算和分配其沸腾质量，计算界面移动距离****/
		/**第1个控制体为气或者气-液**/
		if (sort[1]==2||sort[1]==4)
		{
			n_liquid=1;
			while (n_liquid<=n_liquid_total)
			{
				/*液塞左界面在蒸发段1*/
				if (n_trans[2*n_liquid-1]>=n_heat_start1&&n_trans[2*n_liquid-1]<=n_heat_end1)
				{
					/*右界面超过蒸发段1*/
					if (n_trans[2*n_liquid]>n_heat_end1)
					{
						n=n_trans[2*n_liquid-1];	//从液塞左界面到蒸发段1终点
						while (n<=n_heat_end1)	
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//左界面沸腾传质分配给左侧气泡
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面因沸腾移动的距离
						/*右界面在蒸发段2*/
						if ((n_trans[2*n_liquid]>=n_heat_start2)&&(n_trans[2*n_liquid]<=n_heat_end2))
						{
							n=n_heat_start2;	//从蒸发段2起点到右界面左侧控制体
							while (n<=n_trans[2*n_liquid]-1)
							{
								if (fabs(v_lright[n][1])<=v_minset)
								{
									h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
								}
								else
								{
									h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
								}
								if (T_w[n][1]>=T_lright[n][1])
								{
									Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到右界面
									Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
									if (sort[n]==1)
									{
										Tr_vl_left_each[n]=Tr_vl_right_each[n];
									}
								}
								n=n+1;
							}
							n=n_trans[2*n_liquid];	//右界面单独考虑
							if (fabs(v_lleft[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
							}
							if (T_w[n][1]>=T_lleft[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							if (n_liquid==n_liquid_total)//右界面传质分配给右侧气泡
							{
								Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];	
							}
							else
							{
								Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];
							}
							dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面因沸腾传质移动距离
						}
						/*右界面超过蒸发段2*/
						if (n_trans[2*n_liquid]>n_heat_end2)
						{
							n=n_heat_start2;	//蒸发段2起点到终点
							while (n<=n_heat_end2)
							{
								if (fabs(v_lright[n][1])<=v_minset)
								{
									h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
								}
								else
								{
									h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
								}
								if (T_w[n][1]>=T_lright[n][1])
								{
									Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到右界面
									Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
									if (sort[n]==1)
									{
										Tr_vl_left_each[n]=Tr_vl_right_each[n];
									}
								}
								n=n+1;
							}
							if (n_liquid==n_liquid_total)	//右界面沸腾传质分配给右侧气泡
							{
								Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];
							}
							else
							{
								Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];
							}
							dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面因沸腾移动的距离
						}
					}
					/*右界面在蒸发段1*/
					else
					{
						n=n_trans[2*n_liquid-1];	//从左界面到右界面左边控制体
						while (n<=(n_trans[2*n_liquid]-1))
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						n=n_trans[2*n_liquid];	//右界面单独考虑
						if (fabs(v_lleft[n][1])<=v_minset)
						{
							h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
						}
						else
						{
							h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
						}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{
								Tr_vl_right_each[n]=Tr_vl_left_each[n];
							}
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//左界面沸腾传质一半分配给左侧气泡
						if (n_liquid==n_liquid_total)	//左界面沸腾传质一半分配给右侧气泡
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lleft[n_liquid]/2;
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lleft[n_liquid]/2;
						}
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面因沸腾移动的距离
					}	
				}
				/*液塞左界面在蒸发段2*/
				else if (n_trans[2*n_liquid-1]>=n_heat_start2&&n_trans[2*n_liquid-1]<=n_heat_end2)
				{
					/*右界面超过蒸发段2*/
					if (n_trans[2*n_liquid]>n_heat_end2)
					{
						n=n_trans[2*n_liquid-1];	//从左界面到蒸发段2终点
						while (n<=n_heat_end2)
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//左界面沸腾传质分配给左侧气泡
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面因沸腾传质移动的距离
					}
					/*右界面在蒸发段2*/
					else
					{
						n=n_trans[2*n_liquid-1];	//从左界面到右界面左边控制体
						while (n<=(n_trans[2*n_liquid]-1))
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						n=n_trans[2*n_liquid];	//右界面单独考虑
						if (fabs(v_lleft[n][1])<=v_minset)
						{
							h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
						}
						else
						{
							h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
						}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{
								Tr_vl_right_each[n]=Tr_vl_left_each[n];
							}
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//左界面沸腾传质分配一半给左侧气泡
						if (n_liquid==n_liquid_total)	//左界面沸腾传质分配一半给右侧气泡
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lleft[n_liquid]/2;
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lleft[n_liquid]/2;
						}
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面因沸腾传质移动的距离，右界面不移动？？？
					}
				}
				/*液塞左界面在0到蒸发段1起点之间*/
				else if (n_trans[2*n_liquid-1]>=1&&n_trans[2*n_liquid-1]<n_heat_start1)
				{
					/*右界面在蒸发段1*/
					if (n_trans[2*n_liquid]>=n_heat_start1&&n_trans[2*n_liquid]<=n_heat_end1)
					{
						n=n_heat_start1;	//从蒸发段1起点到右界面左边控制体
						while (n<=n_trans[2*n_liquid]-1)
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到右界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;	
						}
						n=n_trans[2*n_liquid];	//右界面单独考虑
						if (fabs(v_lleft[n][1])<=v_minset)
						{
							h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
						}
						else
						{
							h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
						}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{
								Tr_vl_right_each[n]=Tr_vl_left_each[n];
							}
						}
						if (n_liquid==n_liquid_total)	//右界面沸腾传质分配到右侧气泡
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];
						}
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面因沸腾移动的距离
					}
					/*右界面在蒸发段1之后，蒸发段2之前*/
					else if (n_trans[2*n_liquid]>n_heat_end1&&n_trans[2*n_liquid]<n_heat_start2)
					{
						n=n_heat_start1;	//蒸发段1起点到终点
						while (n<=n_heat_end1)
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//左界面沸腾传质分配一般给左侧气泡
						if (n_liquid==n_liquid_total)	//左界面沸腾传质分配一半给右侧气泡
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lleft[n_liquid]/2;
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lleft[n_liquid]/2;
						}
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面因沸腾移动的距离，右界面不移动？？？
					}
					/*右界面在蒸发段2*/
					else if (n_trans[2*n_liquid]>=n_heat_start2&&n_trans[2*n_liquid]<=n_heat_end2)
					{
						n=n_heat_start1;//蒸发段1起点到终点
						while (n<=n_heat_end1)	
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//左界面沸腾传质分配给左侧气泡
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面因沸腾移动的距离
						n=n_heat_start2;	//蒸发段2起点到右界面
						while (n<=n_trans[2*n_liquid])
						{
							if (fabs(v_lleft[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
							}
							if (T_w[n][1]>=T_lleft[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//沸腾传质累加到右界面
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						if (n_liquid==n_liquid_total)	//右界面沸腾传质分配给右侧气泡
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];
						}
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面因沸腾移动的距离
					}
					/*右界面超过蒸发段2*/
					else if (n_trans[2*n_liquid]>n_heat_end2)
					{
						n=n_heat_start1;	//蒸发段1起点到终点
						while (n<=n_heat_end1)
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//左界面沸腾传质分配给左侧气泡
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面因沸腾移动的距离
						n=n_heat_start2;	//蒸发段2起点到终点
						while (n<=n_heat_end2)
						{
							if (fabs(v_lleft[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
							}
							if (T_w[n][1]>=T_lleft[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//沸腾传质累加到右界面
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						if (n_liquid==n_liquid_total)	//右界面沸腾传质分配给右侧气泡
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];
						}
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面因沸腾移动的距离
					}
				}
				/*液塞左界面在蒸发段1之后，蒸发段2之前*/
				else if (n_trans[2*n_liquid-1]>n_heat_end1&&n_trans[2*n_liquid-1]<n_heat_start2)
				{
					/*右界面在蒸发段2*/
					if (n_trans[2*n_liquid]>=n_heat_start2&&n_trans[2*n_liquid]<=n_heat_end2)
					{
						n=n_heat_start2;	//从蒸发段2起点到右界面
						while (n<=n_trans[2*n_liquid])
						{
							if (fabs(v_lleft[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
							}
							if (T_w[n][1]>=T_lleft[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//沸腾传质累加到右界面
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						if (n_liquid==n_liquid_total)	//右界面沸腾传质分配给右侧气泡
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];
						}
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面因沸腾移动的距离
					}
					/*右界面超过蒸发段2*/
					else if (n_trans[2*n_liquid]>n_heat_end2)
					{
						n=n_heat_start2;	//蒸发段2起点到终点
						while (n<=n_heat_end2)
						{
							if (fabs(v_lleft[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
							}
							if (T_w[n][1]>=T_lleft[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//沸腾传质累加到右界面
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid]/2;	//右界面沸腾传质分配一半给左侧气泡
						if (n_liquid==n_liquid_total)	//右界面沸腾传质分配一半给右侧气泡
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid]/2;
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid]/2;
						}
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面因沸腾移动的距离
					}
				}
				n_liquid=n_liquid+1;
			}
		}
		
		/**第1个控制体为液或者液-气或液-气-液**/
		else
		{
			/**第1个液塞跨过坐标原点，须单独考虑，其左右界面的沸腾和移动分开考虑**/
			n_liquid=1; 
			/*第1个液塞右界面在蒸发段1*/
			if (n_trans[1]>=n_heat_start1&&n_trans[1]<=n_heat_end1)
			{
				n=n_heat_start1;	//蒸发段1起点到右界面
				while (n<=n_trans[1])
				{
					if (fabs(v_lleft[n][1])<=v_minset)
					{
						h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
					}
					else
					{
						h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
					}
					if (T_w[n][1]>=T_lleft[n][1])
					{
						Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//沸腾传质累加到右界面
						Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
						if (sort[n]==1)
						{
							Tr_vl_right_each[n]=Tr_vl_left_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//右界面沸腾传质分配给右侧气泡
				dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面因沸腾移动的距离
			}
			/*第1个液塞右界面在蒸发段1之后，蒸发段2之前*/
			else if (n_trans[1]>n_heat_end1&&n_trans[1]<n_heat_start2)
			{
				n=n_heat_start1;	//蒸发段1起点到终点
				while (n<=n_heat_end1)
				{
					if (fabs(v_lleft[n][1])<=v_minset)
					{
						h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
					}
					else
					{
						h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
					}
					if (T_w[n][1]>=T_lleft[n][1])
					{
						Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//沸腾传质累加到右界面
						Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
						if (sort[n]==1)
						{
							Tr_vl_right_each[n]=Tr_vl_left_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid]/2;	//右界面沸腾传质分配一半给右侧气泡
				Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lright[n_liquid]/2;	//右界面沸腾传质分配一半给左侧气泡
				dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面因沸腾移动的距离
			}
			/*第1个液塞右界面在蒸发段2*/
			else if (n_trans[1]>=n_heat_start2&&n_trans[1]<=n_heat_end2)
			{
				n=n_heat_start1;	//蒸发段1起点到终点
				while (n<=n_heat_end1)
				{
					if (fabs(v_lleft[n][1])<=v_minset)
					{
						h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
					}
					else
					{
						h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
					}
					if (T_w[n][1]>=T_lleft[n][1])
					{
						Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//沸腾传质累加到左界面
						Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
						if (sort[n]==1)
						{
							Tr_vl_right_each[n]=Tr_vl_left_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid];	//左界面沸腾传质分配给左侧气泡
				dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面因沸腾移动的距离
				n=n_heat_start2;	//蒸发段2起点到右界面
				while (n<=n_trans[1])
				{
					if (fabs(v_lleft[n][1])<=v_minset)
					{
						h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
					}
					else
					{
						h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
					}
					if (T_w[n][1]>=T_lleft[n][1])
					{
						Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//沸腾传质累加到右界面
						Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
						if (sort[n]==1)
						{
							Tr_vl_right_each[n]=Tr_vl_left_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//右界面沸腾传质分配给右侧气泡
				dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面因沸腾移动的距离
			}
			
			/*第1个液塞左界面在蒸发段1*/
			if (n_trans[2*n_liquid_total]>=n_heat_start1&&n_trans[2*n_liquid_total]<=n_heat_end1)
			{
				n=n_trans[2*n_liquid_total];	//左界面到蒸发段1终点
				while (n<=n_heat_end1)
				{
					if (fabs(v_lright[n][1])<=v_minset)
					{
						h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
					}
					else
					{
						h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
					}
					if (T_w[n][1]>=T_lright[n][1])
					{
						Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
						Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
						if (sort[n]==1)
						{
							Tr_vl_left_each[n]=Tr_vl_right_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid];	//左界面沸腾传质分配给左侧气泡
				dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面因沸腾移动的距离
				n=n_heat_start2;	//蒸发段1起点到终点
				while (n<=n_heat_end2)
				{
					if (fabs(v_lright[n][1])<=v_minset)
					{
						h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
					}
					else
					{
						h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
					}
					if (T_w[n][1]>=T_lright[n][1])
					{
						Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到右界面
						Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
						if (sort[n]==1)
						{
							Tr_vl_left_each[n]=Tr_vl_right_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//右界面沸腾传质分配给右侧气泡
				dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
			}
			/*第1个液塞左界面在蒸发段1之后，蒸发段2之前*/
			else if (n_trans[2*n_liquid_total]>n_heat_end1&&n_trans[2*n_liquid_total]<n_heat_start2)
			{
				n=n_heat_start2;	//蒸发段2起点到终点
				while (n<=n_heat_end2)
				{
					if (fabs(v_lright[n][1])<=v_minset)
					{
						h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
					}
					else
					{
						h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
					}
					if (T_w[n][1]>=T_lright[n][1])
					{
						Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
						Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
						if (sort[n]==1)
						{
							Tr_vl_left_each[n]=Tr_vl_right_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid]/2;	//左界面沸腾传质分配一半给左侧气泡
				Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//左界面沸腾传质分配一半给右侧气泡
				dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离，右界面不移动？？？
			}
			/*第1个液塞左界面在蒸发段2*/
			else if (n_trans[2*n_liquid_total]>=n_heat_start2&&n_trans[2*n_liquid_total]<=n_heat_end2)
			{
				n=n_trans[2*n_liquid_total];	//左界面到蒸发段2终点
				while (n<=n_heat_end2)
				{
					if (fabs(v_lright[n][1])<=v_minset)
					{
						h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
					}
					else
					{
						h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
					}
					if (T_w[n][1]>=T_lright[n][1])
					{
						Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
						Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
						if (sort[n]==1)
						{
							Tr_vl_left_each[n]=Tr_vl_right_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid];	//左界面沸腾传质分配给左侧气泡
				dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
			}
			
			/**第2个及其后液塞**/
			n_liquid=2; 
			while (n_liquid<=n_liquid_total)
			{
				/*左界面在蒸发段1*/
				if (n_trans[2*n_liquid-2]>=n_heat_start1&&n_trans[2*n_liquid-2]<=n_heat_end1)
				{
					/*右界面超过蒸发段1*/
					if (n_trans[2*n_liquid-1]>n_heat_end1)
					{
						n=n_trans[2*n_liquid-2];	//左界面到蒸发段1终点
						while (n<=n_heat_end1)
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//左界面沸腾传质分配给左侧气泡
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
						/*右界面在蒸发段2*/
						if ((n_trans[2*n_liquid-1]>=n_heat_start2)&&(n_trans[2*n_liquid-1]<=n_heat_end2))
						{
							n=n_heat_start2;	//蒸发段2起点到右界面左边控制体
							while (n<=n_trans[2*n_liquid-1]-1)
							{
								if (fabs(v_lright[n][1])<=v_minset)
								{
									h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
								}
								else
								{
									h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
								}
								if (T_w[n][1]>=T_lright[n][1])
								{
									Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到右界面
									Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
									if (sort[n]==1)
									{
										Tr_vl_left_each[n]=Tr_vl_right_each[n];
									}
								}
								n=n+1;
							}
							n=n_trans[2*n_liquid-1];	//右界面单独考虑？？？
							if (fabs(v_lleft[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
							}
							if (T_w[n][1]>=T_lleft[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//右界面沸腾传质分配到右侧气泡
							dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
						}
						/*右界面超过蒸发段2*/
						if (n_trans[2*n_liquid-1]>n_heat_end2)
						{
							n=n_heat_start2;	//蒸发段2起点到终点
							while (n<=n_heat_end2)
							{
								if (fabs(v_lright[n][1])<=v_minset)
								{
									h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
								}
								else
								{
									h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
								}
								if (T_w[n][1]>=T_lright[n][1])
								{
									Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到右界面
									Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
									if (sort[n]==1)
									{
										Tr_vl_left_each[n]=Tr_vl_right_each[n];
									}
								}
								n=n+1;
							}
							Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//右界面沸腾传质分配给右侧气泡
							dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
						}
					}
					/*右界面在蒸发段1*/
					else
					{
						n=n_trans[2*n_liquid-2];	//左界面到右界面左边控制体
						while (n<=(n_trans[2*n_liquid-1]-1))
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						n=n_trans[2*n_liquid-1];	//右界面单独考虑
						if (fabs(v_lleft[n][1])<=v_minset)
						{
							h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
						}
						else
						{
							h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
						}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{
								Tr_vl_right_each[n]=Tr_vl_left_each[n];
							}
						}
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid]/2;	//左界面沸腾传质分配一半给左侧气泡
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//左界面沸腾传质分配一半给右侧气泡
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离，右界面不动？？？
					}	
				}
				/*左界面在蒸发段2*/
				else if (n_trans[2*n_liquid-2]>=n_heat_start2&&n_trans[2*n_liquid-2]<=n_heat_end2)
				{
					/*右界面超过蒸发段2*/
					if (n_trans[2*n_liquid-1]>n_heat_end2)
					{
						n=n_trans[2*n_liquid-2];	//左界面到蒸发段2终点
						while (n<=n_heat_end2)
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//左界面沸腾传质分配给左侧气泡
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
					}
					/*右界面在蒸发段2*/
					else
					{
						n=n_trans[2*n_liquid-2];	//左界面到右界面左边控制体
						while (n<=(n_trans[2*n_liquid-1]-1))
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						n=n_trans[2*n_liquid-1];	//右界面单独考虑
						if (fabs(v_lleft[n][1])<=v_minset)
						{
							h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
						}
						else
						{
							h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
						}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{
								Tr_vl_right_each[n]=Tr_vl_left_each[n];
							}
						}
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid]/2;	//左界面沸腾传质分配一半给左侧气泡
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//左界面沸腾传质分配一半给右侧气泡
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离，右界面不移动？？？
					}
				}
				/*左界面在1和蒸发段1之间*/
				else if (n_trans[2*n_liquid-2]>=1&&n_trans[2*n_liquid-2]<n_heat_start1)
				{
					/*右界面在蒸发段1*/
					if (n_trans[2*n_liquid-1]>=n_heat_start1&&n_trans[2*n_liquid-1]<=n_heat_end1)
					{
						n=n_heat_start1;	//蒸发段1起点到右界面左边控制体
						while (n<=n_trans[2*n_liquid-1]-1)
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到右界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;	
						}
						n=n_trans[2*n_liquid-1];	//右界面单独考虑
						if (fabs(v_lleft[n][1])<=v_minset)
						{
							h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
						}
						else
						{
							h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
						}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
							if (sort[n]==1)
							{
								Tr_vl_right_each[n]=Tr_vl_left_each[n];
							}
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//右界面沸腾传质分配给右侧气泡
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
					}
					/*右界面在蒸发段1之后，蒸发段2之前*/
					else if (n_trans[2*n_liquid-1]>n_heat_end1&&n_trans[2*n_liquid-1]<n_heat_start2)
					{
						n=n_heat_start1;	//蒸发段1起点到终点
						while (n<=n_heat_end1)
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid]/2;	//左界面沸腾传质分配一半给左侧气泡
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//左界面沸腾传质分配一半给右侧气泡
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);		//左界面移动距离，右界面不移动？？？
					}
					/*右界面在蒸发段2*/
					else if (n_trans[2*n_liquid-1]>=n_heat_start2&&n_trans[2*n_liquid-1]<=n_heat_end2)
					{
						n=n_heat_start1;	//蒸发段1起点到终点
						while (n<=n_heat_end1)
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//左界面沸腾传质分配给左侧气泡
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
						n=n_heat_start2;	//蒸发段2起点到右界面
						while (n<=n_trans[2*n_liquid-1])
						{
							if (fabs(v_lleft[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
							}
							if (T_w[n][1]>=T_lleft[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//沸腾传质累加到右界面
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//右界面沸腾传质分配给右侧气泡
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
					}
					/*右界面超过蒸发段2*/
					else if (n_trans[2*n_liquid-1]>n_heat_end2)
					{
						n=n_heat_start1;	//蒸发段1起点到终点
						while (n<=n_heat_end1)
						{
							if (fabs(v_lright[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
							}
							if (T_w[n][1]>=T_lright[n][1])
							{
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//沸腾传质累加到左界面
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//左界面沸腾传质分配给左侧气泡
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//左界面移动距离
						n=n_heat_start2;	//蒸发段2起点到终点
						while (n<=n_heat_end2)
						{
							if (fabs(v_lleft[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
							}
							if (T_w[n][1]>=T_lleft[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//沸腾传质累加到右界面
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//右界面沸腾传质分配给右侧气泡
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面动距离
					}
				}
				/*左界面在蒸发段1之后，蒸发段2之前*/
				else if (n_trans[2*n_liquid-2]>n_heat_end1&&n_trans[2*n_liquid-2]<n_heat_start2)
				{
					/*右界面在蒸发段2*/
					if (n_trans[2*n_liquid-1]>=n_heat_start2&&n_trans[2*n_liquid-1]<=n_heat_end2)
					{
						n=n_heat_start2;	//蒸发段2起点到右界面
						while (n<=n_trans[2*n_liquid-1])
						{
							if (fabs(v_lleft[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
							}
							if (T_w[n][1]>=T_lleft[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//沸腾传质累加到右界面
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//右界面沸腾传质分配给右侧气泡
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离
					}
					/*右界面超过蒸发段2*/
					else if (n_trans[2*n_liquid-1]>n_heat_end2)
					{
						n=n_heat_start2;	//蒸发段2起点到终点
						while (n<=n_heat_end2)
						{
							if (fabs(v_lleft[n][1])<=v_minset)
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);
							}
							else
							{
								h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);
							}
							if (T_w[n][1]>=T_lleft[n][1])
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//沸腾传质累加到右界面
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lright[n_liquid]/2;	//右界面沸腾传质分配一半给左侧气泡
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid]/2;	//右界面沸腾传质分配一半给右侧气泡
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//右界面移动距离，左界面不移动？？？
					}
				}
				n_liquid=n_liquid+1;
			}
		}
		
		/*重新设定相变标记*/
		for (n=1;n<=length;n=n+1)
		{
			Tr_vl_lright_sig[n]=0;
			Tr_vl_lleft_sig[n]=0;
			if (heat_sig[n]==1)
			{
				if (T_lright[n][1]!=0&&(T_w[n][1]>=T_lright[n][1]))
				{
					Tr_vl_lright_sig[n]=1;
				}

				if (T_lleft[n][1]!=0&&(T_w[n][1]>=T_lleft[n][1]))
				{
					Tr_vl_lleft_sig[n]=1;
				}
			}
		}
		/*计算液塞总长度*/
		double l_liquid_total=0;	//液体总长度
		n=1;
		while (n<=length)
		{
			if (sort[n]==1)
			{
				l_liquid_total=l_liquid_total+x_lleft[n][1];
			}
			else
			{
				l_liquid_total=l_liquid_total+x_lleft[n][1]+x_lright[n][1];
			}
			n=n+1;
		}
		/*检查控制体内液体长度是否出现‘无效数字’或‘无穷大’*/
		int special_nan=0;	//出现无效数字或无穷大的标记
		for (n=1;n<=length;n=n+1)
		{
			if (_isnan(x_lleft[n][1])!=0||_finite(x_lleft[n][1])==0||_isnan(x_lright[n][1])!=0||_finite(x_lright[n][1])==0)
			{
				special_nan=1;
				break;	//中止循环for (n=1;n<=length;n=n+1)
			}
		}
		/*检查控制体内壁温是否出现‘无效数字’或‘无穷大’*/
		for (n=1;n<=length;n=n+1)
		{
			if (_isnan(T_w[n][1])!=0||_finite(T_w[n][1])==0)
			{
				special_nan=2;
				break;	//中止循环for (n=1;n<=length;n=n+1
			}
		}
		if (special_nan!=0)
		{
			fid0=fopen("state.txt","a+");	//按照a+模式打开一个名叫state.txt的文件，返回state.txt的文件指针给fid0
			fprintf(fid0,"结束,%d,n=%d\n",special_nan,n); //将"结束，special_nan和n"写入由fid0指出的文件
		}


		/*%%%%%%%输出端%%%%%%%%*/
		if (i%print_fre==0 )	//时间步长数量满足输出要求
		{
			kk=kk+1;   
						
			fid0=fopen("state.txt","a+");	//按照a+模式打开一个名叫state.txt的文件，返回state.txt的文件指针给fid0
			fprintf(fid0,"第%d组数据计算完毕\n",kk);	//将"第kk组数据计算完毕"写入由fid0指出的文件
			fclose(fid0);	//关闭由fid0指出的文件,返回操作结果，0或EOF
			printf("第%d组数据计算完毕\n",kk);	//在控制台输出"第kk组数据计算完毕"

			fid1=fopen("ro_v.txt","a+"); //按照a + 模式打开一个名叫ro_v.txt的文件，返回ro_v.txt的文件指针给fid1
			fid2=fopen("T_w.txt","a+");	//按照a + 模式打开一个名叫T_w.txt的文件，返回T_w.txt的文件指针给fid2
			fid3=fopen("trans_point.txt","a+");
			fid4=fopen("P_v.txt","a+");
			fid5=fopen("T_e.txt","a+");
			fid6=fopen("x_l.txt","a+");
			fid7=fopen("T_e2.txt","a+");
			fid8=fopen("T_v.txt","a+");
			fid22=fopen("T_sat.txt","a+");
			fid9=fopen("v_l.txt","a+");
			fid10=fopen("T_l.txt","a+");
			
			fid13=fopen("T_vb.txt","a+");
			fid23 = fopen("T_sat_vb.txt", "a+");
			fid14=fopen("P_vb.txt","a+");
			fid15=fopen("ro_vb.txt","a+");
			fid16=fopen("P_l.txt","a+");
			fid18=fopen("bubble_gen_total.txt","a+");
			fid19=fopen("P_place.txt","a+");
			fid20=fopen("h_b.txt","a+");
			fid21=fopen("h_l.txt","a+");
			
			
			fprintf(fid1,"ro_v[n]数据组数是：%d ",kk);	//将"数据组数是：kk"写入由fid1指出的文件ro_v.txt
			fprintf(fid2,"T_w[n]数据组数是：%d ",kk);
			fprintf(fid3,"trans_point数据组数是：%d ",kk);
			fprintf(fid4,"P_v[n]数据组数是：%d ",kk);
			fprintf(fid5,"Te数据组数是：%d ",kk);
			fprintf(fid6,"x_llert[n],x_lright[n]数据组数是：%d ",kk);
			fprintf(fid7, "Te2数据组数是：%d ", kk);
			fprintf(fid8,"T_v[n]数据组数是：%d ",kk);
			fprintf(fid22, "T_sat[n]数据组数是：%d ", kk);
			fprintf(fid9,"v_lleft[n],v_lright[n]数据组数是：%d ",kk);
			fprintf(fid10,"T_lleft[n],T_lright[n]数据组数是：%d ",kk);
			
			fprintf(fid13,"T_vb[n_bubble]数据组数是：%d ",kk);
			fprintf(fid23, "T_sat_vb[n_bubble]数据组数是：%d ", kk);
			fprintf(fid14,"P_vb[n_bubble]数据组数是：%d ",kk);
			fprintf(fid15,"ro_vb[n_bubble]数据组数是：%d ",kk);
			fprintf(fid16,"P_l[n]数据组数是：%d ",kk);
			fprintf(fid18,"数据组数是：%d ",kk);
			fprintf(fid20,"h_b[n]数据组数是：%d ",kk);
			fprintf(fid21,"h_lleft[n],h_lright[n]数据组数是：%d ",kk);
			
			n_bubble=1;	//输出每个气泡的温度、压力、密度
			while (n_bubble<=n_bubble_total)
			{
				fprintf(fid13,"T_vb[%d]=%f,",n_bubble,T_vb[n_bubble][1]);	//将"T_vb[n_bubble]=T_vb[n_bubble][1],"写入由fid13指出的文件T_vb.txt
				fprintf(fid23, "T_sat_vb[%d]=%f,", n_bubble, T_sat_vb[n_bubble][1]);	//将"T_sat_vb[n_bubble]=T_sat_vb[n_bubble][1],"写入由fid23指出的文件T_sat_vb.txt
				fprintf(fid14,"P_vb[%d]=%f,",n_bubble,P_vb[n_bubble][1]);	//将"P_vb[n_bubble]=P_vb[n_bubble][1],"写入由fid14指出的文件P_vb.txt
				fprintf(fid15,"ro_vb[%d]=%f,",n_bubble,ro_vb[n_bubble][1]);	//将"ro_vb[n_bubble]=ro_vb[n_bubble][1],"写入由fid15指出的文件ro_vb.txt
				n_bubble=n_bubble+1;
			}
			fprintf(fid13,"\n");	//将换行写入由fid13指出的文件T_vb.txt
			fprintf(fid23, "\n");	//将换行写入由fid23指出的文件T_sat_vb.txt
			fprintf(fid14,"\n");
			fprintf(fid15,"\n");

			
			location = 1;       /*%数据坐标，相当于n*/
			while (location<=length)	//输出每个控制体的参数，写入txt文件
			{
				fprintf(fid1,"%f ",ro_v[location][1]);	//%f后面的空格用于导入Excel时分列
				fprintf(fid2,"%f ",T_w[location][1]);
				fprintf(fid4,"%f ",P_v[location][1]);
				fprintf(fid6,"%f,%f ",x_lleft[location][1],x_lright[location][1]);
				fprintf(fid8,"%f ",T_v[location][1]);
				fprintf(fid22, "%f ", T_sat[location][1]);
				fprintf(fid9,"%f,%f ",v_lleft[location][1],v_lright[location][1]);
				fprintf(fid10,"%f,%f ",T_lleft[location][1],T_lright[location][1]);
				fprintf(fid16,"%f ",P_l[location][2]);
				fprintf(fid20,"%f ",h_b[location]);
				fprintf(fid21,"%f,%f ",h_lleft[location],h_lright[location]);
				
				/*写入蒸发段1和2的壁温，以及绝热段中部的压力*/
				if (location==155)
				{
					fprintf(fid5,"T_w[155]=%f\n",T_w[155][1]);
				}
				if (location==575)
				{
					fprintf(fid7,"T_w[575]=%f\n",T_w[575][1]);
				}
				if (location==260)
				{
					if (sort[location]==1)
					{
						fprintf(fid19,"sort[260]=1,P[260]=%f\n",P_l[260][2]);
					}
					else
					{
						fprintf(fid19,"sort[260]=%d,P[260]=%f\n",sort[260],P_v[260][2]);
					}
				}

				/*写入气液界面所在的控制体坐标，及界面控制体的液体长度*/
				if (sort[location]==3)
				{
					fprintf(fid3,"n=%d sort[%d]=3,x_lleft=%f,x_lright=%f ", location, location,x_lleft[location][1],x_lright[location][1]);
				}

				if (sort[location]==4)
				{
					fprintf(fid3,"n=%d sort[%d]=4,x_lleft=%f,x_lright=%f ", location, location,x_lleft[location][1],x_lright[location][1]);
				}
				if (sort[location]==8)
				{
					fprintf(fid3,"n=%d sort[%d]=8,x_lleft=%f,x_lright=%f ", location, location,x_lleft[location][1],x_lright[location][1]);
				}
				location=location+1;
			}

			/*遍历n再写入换行*/
			fprintf(fid1, "\n");
			fprintf(fid2, "\n");
			fprintf(fid3, "\n");
			fprintf(fid4, "\n");
			fprintf(fid6, "\n");
			fprintf(fid8, "\n");
			fprintf(fid22, "\n");
			fprintf(fid9, "\n");
			fprintf(fid10, "\n");
			fprintf(fid16, "\n");
			fprintf(fid20, "\n");
			fprintf(fid21, "\n");

			fprintf(fid18,"bubble_gen_total=%d\n",bubble_gen_total);
		
			fclose(fid1);
			fclose(fid2);
			fclose(fid3);
			fclose(fid4);
			fclose(fid5);
			fclose(fid6);
			fclose(fid7);
			fclose(fid8);
			fclose(fid9);
			fclose(fid10);
			fclose(fid13);
			fclose(fid14);
			fclose(fid15);
			fclose(fid16);
			fclose(fid18);
			fclose(fid19);
			fclose(fid20);
			fclose(fid21);
			fclose(fid22);
			fclose(fid23);

		}
		/*写入显热，潜热，液体长度*/
		if ((i%heat_print_fre)==0)	//时间步长数量满足输出要求
		{
			sensible_c_total=sensible_c_total/heat_print_fre;	//前面累加，此处对时间步长数做平均
			sensible_e_total=sensible_e_total/heat_print_fre;
			latent_c_total=latent_c_total/heat_print_fre;
			latent_e_total=latent_e_total/heat_print_fre;
			fid11=fopen("sensible.txt","a+");
			fid12=fopen("latent.txt","a+");
			fid17=fopen("liquid_total.txt","a+");
			fprintf(fid11,"时间(s)：%f,sensible_c_total=%f,sensible_e_total=%f\n",i/heat_print_fre*0.001,sensible_c_total,sensible_e_total);
			fprintf(fid12,"时间(s)：%f,latent_c_total=%f,latent_e_total=%f\n",i/heat_print_fre*0.001,latent_c_total,latent_e_total);
			fprintf(fid17,"时间(s)：%f,liquid_total=%f\n",i/heat_print_fre*0.001,l_liquid_total);
			fclose(fid11);
			fclose(fid12);
			fclose(fid17);
			sensible_c_total=0;	//重置
			sensible_e_total=0;
			latent_c_total=0;
			latent_e_total=0;
		}
		else  //时间步长不满足输出要求
		{
			sensible_c_total=sensible_c_total+sensible_c;	//将这一时刻的冷凝段显热累加
			sensible_e_total=sensible_e_total+sensible_e;
			latent_c_total=latent_c_total+latent_c;
			latent_e_total=latent_e_total+latent_e;
		}

		n=1;
		while (n<=length)
		{
			P_l[n][2]=0;	//在主循环内重置液体压力
			n=n+1;
		}

		i=i+1;	//下一个时间节点
	}
	

	fid0=fopen("state.txt","a+");	//按照a+模式打开一个名叫state.txt的文件，返回state.txt的文件指针给fid0
	if (codenumber==1||codenumber==2)
		{fprintf(fid0,"Out of range, Codenumber=%d, n=%d, i=%d",codenumber,n,i);}	//将codenumber,n,i写入由fid0指出的文件
	else
		{fprintf(fid0,"The end\n");}	//将“The end\n”写入由fid0指出的文件
	fclose(fid0);	//关闭由fid0指出的文件,返回操作结果，0或EOF

	printf("The end");	//在屏幕输出The end
	getchar();	//从计算机终端输入一个字符。作用是程序结束时防止命令窗消失，不能观看运行结果，设置一个等待
	return 0;
}


/*****************************************************************
 * obtain liquid dynamic viscosity at a given temperature, Pa-s  *
 * Liquid hydrogen properties fitting.EES                        *
 *****************************************************************/
double mu_l_given_T(double T)
{
	double output = 0.0000882776 - 0.00000729164*T + 2.33925E-07*pow(T, 2) - 2.72711E-09*pow(T, 3);
	return output;
}

/***************************************************
 * obtain cp of LH2 at a given temperature, J/kg-K *
 * Liquid hydrogen properties fitting.EES          *
 ***************************************************/
double cp_l_given_T(double T)
{
	double output = -2.92658E+06 + 646211 * T - 56670.1*pow(T, 2) + 2473.7*pow(T, 3) - 53.7465*pow(T, 4) + 0.465581*pow(T, 5);
	return output;
}

/***************************************************
 * obtain cv of LH2 at a given temperature, J/kg-K *
 * Liquid hydrogen properties fitting.EES          *
 ***************************************************/
double cv_l_given_T(double T)
{
	double output = -38270.6 + 9394.37*T - 815.337*pow(T, 2) + 35.5592*pow(T, 3) - 0.774024*pow(T, 4) + 0.00672672*pow(T, 5);
	return output;
}

/********************************************************************
 * obtain thermal conductivity of LH2 at a given temperature, W/m-K *
 * Liquid hydrogen properties fitting.EES                           *
 ********************************************************************/
double lamt_l_given_T(double T)
{
	double output = 0.0259393 + 0.00833281*T - 0.000267213*pow(T, 2) + 0.00000229541*pow(T, 3);
	return output;
}

/**********************************************************
 * obtain Prandtl number of LH2 at a given temperature, - *
 * Liquid hydrogen properties fitting.EES                 *
 **********************************************************/
double Pr_l_given_T(double T)
{
	double output = -110.962 + 25.3724*T - 2.26206*pow(T, 2) + 0.099788*pow(T, 3) - 0.00218533*pow(T, 4) + 0.0000190604*pow(T, 5);
	return output;
}

/*******************************************************
 * obtain density of LH2 at a given temperature, kg/m3 *
 * Liquid hydrogen properties fitting.EES              *
 *******************************************************/
double ro_l_given_T(double T)
{
	double output = 128.792 - 6.49526*T + 0.274511*pow(T, 2) - 0.00468565*pow(T, 3);
	return output;
}

 /******************************************************************
  * obtain P_sat at a given temperature, fitted from RKS data, Pa  *
  * Psat from Redilich-Kwong-Soave(RKS)EOS.EES                     *
  ******************************************************************/
double p_sat_given_T(double T)
{
	double output = -19613 + 27548.8*T - 3374.38*pow(T, 2) + 112.977*pow(T, 3);
	return output;
}

/****************************************************************
 * obtain T_sat at a given temperature, fitted from RKS data, K *
 * Psat from Redilich-Kwong-Soave(RKS)EOS.EES                   *
 ****************************************************************/
double T_sat_given_p(double p)
{
	double output = 16.0673 + 0.000056766*p - 1.58057E-10*pow(p, 2) + 2.82160E-16*pow(p, 3) - 2.54034E-22*pow(p, 4) + 8.87606E-29*pow(p, 5);
	return output;
}

/**********************************************************************
 * obtain rho_sat at a given temperature, fitted from RKS data, kg/m3 *
 * Psat from Redilich-Kwong-Soave(RKS)EOS.EES                         *
 **********************************************************************/
double ro_sat_given_T(double T)
{
	double output = 96.8941 - 17.4254*T + 1.17089*pow(T, 2) - 0.0351513*pow(T, 3) + 0.000409864*pow(T, 4);
	return output;
}

/***********************************************************************************
 * obtain vaporization enthalpy at a given temperature, fitted from RKS data, J/kg *
 * Psat from Redilich-Kwong-Soave(RKS)EOS.EES                                      *
 ***********************************************************************************/
double h_fg_given_T(double T)
{
	double output = -5.64364E+06 + 1.01044E+06*T - 61601.1*pow(T, 2) + 1662.36*pow(T, 3) - 17.0694*pow(T, 4);
	return output;
}

/********************************************************************************
 * obtain internal energy from specivic volume and temperature, using RKS, J/kg *
 * reference state: 16.47m3/kg, 40K                                             *
 * details in "u from RKS-analytical.EES"                                       *
 ********************************************************************************/
double u_given_vt(double volume, double mass, double temp)
{
	double v = volume / mass;	//specific volume, m3/kg
	double Tr = temp / 33.145f;	//reduced temperature, -
	double m = 0.127348;	//paremeter used in RKS
	double alpha = pow(1 + m * (1 - sqrt(Tr)), 2);	//paremeters used in RKS
	double c = 6162.24602;	//that is a/alpha in RKS
	double b = 0.009136094;	//parameer used in RKS
	double u_ref = 505138.324;	//reference internal energy, J/kg
	double v_ref = 16.47;	//reference specific volume, m3/kg
	double integral_T = -251136.681 + 6877.58*temp - 110.209 / 2.0*pow(temp, 2) + 6.99154 / 3.0*pow(temp, 3) - 0.219446 / 4.0*pow(temp, 4) + 0.00340588 / 5.0*pow(temp, 5) - 0.0000208872 / 6.0*pow(temp, 6);
	double integral_v = c * (m*sqrt(alpha)*sqrt(Tr) + alpha) / b * (log(v / (v + b)) - log(v_ref / (v_ref + b)));
	double u = u_ref + integral_T + integral_v;
	return u;
}

/*****************************************************************************
 *obtain temperature from internal energy and specific volume, K             *
 *using the previous temperature T1 as the guess value of new temperature T2 *
 *****************************************************************************/
double T_given_uv(double u, double volume, double mass, double T1)
{
	int MAX_ITERS = 1000;      //maximun iteration times
	double MAX_DELTA = 0.0001;  //maximun internal energy difference
	double step = 0.01;       //step to guess temperature
	double guess = T1;
	double mark = u;
	double delta = fabs(u_given_vt(volume, mass, guess) - mark);
	int i = 0;
	for (i = 0; (i < MAX_ITERS) && (delta > MAX_DELTA); i++)
	{
		double delta_a = fabs(u_given_vt(volume, mass, guess + step) - mark);
		double delta_b = fabs(u_given_vt(volume, mass, guess - step / 2.0) - mark);
		if (delta_a < delta_b)
		{
			delta = delta_a;
		}
		else
		{
			delta = delta_b;
			step = -step / 2.0;
		}
		guess += step;
	}
	return guess;
}
