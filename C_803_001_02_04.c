#include <stdio.h>
#include <math.h>
#include <float.h>


#define pi 3.141592653589
#define length 840	//�ܳ���840mm�����ȶ�100mm��������������ÿ��110mm

/*����ԭ��*/
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
	
	/* ���岽�� */
	double dx=1e-3;      /*���벽����m */
	double dt=2e-7;      /*ʱ�䲽��,s*/

	/* ������ʱ��-���������� */
	double T_total=20;	//��ʱ����s
	long nn=T_total/dt;	//��ʱ��ڵ���
	int print_fre=5000;	//����Ĳ��������ÿ5000��ʱ�䲽�����һ�Σ�ÿ0.001s���һ�Σ�
	long heat_print_fre=0.001/dt;     /*����������Ĳ����������0.001s���һ��ƽ����������*/

	/* ���β��� */
	const double d = 2.3e-3;	//PHP�ھ���m
	const double d_out = 3.3e-3;	//PHP�⾶��m
	const double film = 5.2e-5;	//ҺĤ��ȣ�m
	const double A_v = 0.25*pi*pow((d - 2 * film), 2); //���ݽ������m2
	const double A_l=4.15476e-6;	//Һ���������m2���ھ�2.3mm��Բ
	const double A_w = 4.39823e-6;	//�ܱڽ������m2���ھ�2.3mm�⾶3.3mm��Բ��

	const double g=9.8;	//�������ٶȣ�m/s2
	
	/* ���Բ��� */
	const double ro_l = 72.26;	//Һ��19K�ܶȣ�kg/m3
	const double mm = 4124;    //���������峣��EES��J/kg-K
	const double lamt_w = 2.169;                           /*20K,SS304�ȵ��ʣ�W/m-K*/
	const double c_w = 13.45;                            /*20K,SS304����,J/kg-K*/
	const double ro_w = 8072;                           /*20K,SS304�ܱڵ��ܶȣ�kg/m3*/

	const double T_min=18;	//�¶����ޣ�K
	const double T_max=31;	//�¶����ޣ�K
	/* �����¶ȣ�Ԥ������*/
	const double T_out = 20;
	const double T_c = 19.065;	//�̶��������¶ȣ�K
	
	/* ״ֵ̬ */
	int codenumber = 0;	//0��ʾ������1��ʾ��n=1���忪ʼ��Һ���ƶ������ʱ�䲽������2��ʾ��n=1Һ�忪ʼ��Һ���ƶ������ʱ�䲻������
	long bubble_gen_total = 0;	//�²������ݵ�������


    /* ****************����������������*******************/
	int sort[length+1];	//������������飬�����±�Ϊ0-length�������±곤��Ϊlength+1��=1��ʾҺ��=2��ʾ����=3��ʾҺ-����=4��ʾ��-Һ��=8��ʾҺ-��-Һ��=9��ʾҺ��Ҫ�ϲ���=6��ʾ�����ݲ���
	int c_gravity[length+1];	//��������ϵ����ֵΪ-1��0��1
	int plain[length+1];	//������ʽ��ǣ�1��ʾ�趨���£�0��ʾ���¿ɱ�
	float heat_flux[length+1];	//�����ܶȣ�W/m
	int heat_sig[length+1];		//���ȱ�ǣ�1��ʾ�����Σ�0��ʾ���ȶΣ�-1��ʾ������
	float h_v[length+1];	//�����е�������������ϵ����W/m2-K
	float h_c[length+1];	//�����������ϵ������Ϊ0
	double h_b[length+1];	//���ڻ���ϵ����W/m2-K
	int Tr_v_sig[length+1];	//ҺĤ����ǣ�1������0����䣬-1����
	double h_lleft[length+1];	//Һ�嵥���������ϵ����W/m2-K
	double h_lright[length+1];	//Һ�嵥���������ϵ����W/m2-K

	double mu_left[length+1];	//ÿ�������������Һ������ճ�ȣ�Pa-s
	double mu_right[length+1];	//ÿ�����������ұ�Һ������ճ�ȣ�Pa-s
	double c_pl_left[length+1];	//ÿ�������������Һ�����ȣ�J/kg-K
	double c_pl_right[length+1];	//ÿ�����������ұ�Һ�����ȣ�J/kg-K
	double lamt_left[length+1];	//ÿ�������������Һ���ȵ��ʣ�W/m-K
	double lamt_right[length+1];	//ÿ�����������ұ�Һ���ȵ��ʣ�W/m-K
	double h_fg_lleft[length+1];	//ÿ�������������Һ������Ǳ�ȣ�J/kg
	double h_fg_lright[length+1];	//ÿ�����������ұ�Һ������Ǳ�ȣ�J/kg
	double h_fg_v[length+1];	//�����������岿�ֵ�����Ǳ�ȣ���Һ�������Ǳ�����
	double C_l[length+1];	//�������Һ������ϵ��

	double x_lleft[length+1][3];	//ÿ�������������Һ�����ȣ�m
	double x_lright[length+1][3];	//ÿ�����������ұ�Һ�����ȣ�m
	double v_lleft[length+1][3];	//ÿ�������������Һ���ٶȣ�m/s
	double v_lright[length+1][3];	//ÿ�����������ұ�Һ���ٶȣ�m/s
	double T_lleft[length+1][3];	//ÿ�������������Һ���¶ȣ�m/s
	double T_lright[length+1][3];	//ÿ�����������ұ�Һ���¶ȣ�m/s
	double P_l[length+1][3];	//Һ��ѹ����Pa
	double x_v[length+1][3];	//���������������ȣ�m
	double ro_v[length+1][3];	//�������������ܶȣ�kg/m3
	double T_v[length+1][3];	//�������������¶ȣ�K
	double T_sat[length + 1][3];	//��������ҺĤ���汥���¶ȣ�K
	double P_v[length+1][3];	//������������ѹ����Pa

	double T_w[length+1][3];	//�����¶ȣ�K

	/* ��ʼ����Һ�ֲ���10������10��Һ������Һ��50%*/
	int n;	//n��ʾ��������ţ�λ�����꣩����1��840��ÿ�������峤1mm
	{
		sort[0] = 0;	//����ĵ�һ��Ԫ���ò�������ֵ0
		for (n = 1; n <= 20; n = n + 1)
		{sort[n] = 2;}	//����������
		sort[21] = 4;	//��-Һ������
		for (n = 22; n <= 62; n = n + 1)
		{sort[n] = 1;}	//Һ��������
		sort[63] = 3;	//Һ-��������

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
	/* �趨��������ϵ��*/
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
	/* �趨���¹̶�ϵ��*/
	{
		for (n = 1; n <= length; n = n + 1)
		{
			plain[n] = 0;	//���¿ɱ�
		}
		for (n = 311; n <= 420; n = n + 1)
		{
			plain[n] = 1;	//���±��趨��������1
		}
		for (n = 731; n <= 840; n = n + 1)
		{
			plain[n] = 1;	//���±��趨��������2
		}
	}
	/* �趨������(W/m)��������1��2��λ��*/
	float heat_flux_set = 0.4 / 0.11;	//��������������0.8Wƽ�����䵽���������ι�0.11m�ĳ�����
	int n_heat_start1 = 101;	//��������ʼ���������
	int n_heat_end1 = 210;	//��������ֹ���������
	int n_heat_start2 = 521;	//��������ʼ���������
	int n_heat_end2 = 630;	//��������ֹ���������
	/*�趨�����ֿ�����������ܶȺʹ��ȱ��*/
	{
		for (n = 1; n <= 100; n = n + 1)	//���ȶ�1
		{	heat_flux[n] = 0;
			heat_sig[n] = 0;		}
		for (n = 101; n <= 210; n = n + 1)	//������1
		{	heat_flux[n] = heat_flux_set;
			heat_sig[n] = 1;		}
		for (n = 211; n <= 310; n = n + 1)	//���ȶ�2
		{	heat_flux[n] = 0;
			heat_sig[n] = 0;		}
		for (n = 311; n <= 420; n = n + 1)	//������1�������ܶ�Ϊ0
		{	heat_flux[n] = 0;
			heat_sig[n] = -1;		}
		for (n = 421; n <= 520; n = n + 1)	//���ȶ�3
		{	heat_flux[n] = 0;
			heat_sig[n] = 0;		}
		for (n = 521; n <= 630; n = n + 1)	//������2
		{	heat_flux[n] = heat_flux_set;
			heat_sig[n] = 1;		}
		for (n = 631; n <= 730; n = n + 1)	//���ȶ�4
		{	heat_flux[n] = 0;
			heat_sig[n] = 0;		}
		for (n = 731; n <= 840; n = n + 1)	//������2
		{	heat_flux[n] = 0;
			heat_sig[n] = -1;		}
	}
	/* �趨���ݹ��ȶȣ�K*/
	float overheat = 0.25;
	/* �趨���ݲ���ʱ����*/
	double generate_fre = 0.05;	//���ݲ�����ʱ������s
	double generate_inter1 = 0;	//��1�����ݲ�����ļ�ʱ��
	double generate_inter2 = 0;	//��2�����ݲ�����ļ�ʱ��
	double generate_inter3 = 0;	//��3�����ݲ�����ļ�ʱ��
	double generate_inter4 = 0;	//��4�����ݲ�����ļ�ʱ��
	double generate_inter5 = 0;	//��5�����ݲ�����ļ�ʱ��
	double generate_inter6 = 0;	//��6�����ݲ�����ļ�ʱ��
	double generate_inter7 = 0;	//��7�����ݲ�����ļ�ʱ��
	double generate_inter8 = 0;	//��8�����ݲ�����ļ�ʱ��
	/* �趨��������ʧ��С����*/
	double l_disappear = 0.2*dx;
	/* �趨��������������������ϵ��*/
	for (n = 1; n <= length; n = n + 1)
	{
		h_v[n] = 0;
		if (heat_sig[n] == 1)
		{h_v[n] = 2028;}	//��������ϵ��1000W/m2-K
		if (heat_sig[n] == -1)
		{h_v[n] = 2028;}	//��������ϵ��1000W/m2-K
	}
	/* �趨������ȴ���Ȳ�����Ԥ����*/
	for (n = 1; n <= length; n = n + 1)
	{h_c[n] = 0;}
	/* �趨��ʼ�����*/
	for (n = 1; n <= length; n = n + 1)
	{Tr_v_sig[n] = 0;	//��ʼ�趨Ϊ�����
	}
	/*���㽻���������λ��*/
	int transfer_point_number;	//������¼��Һ�������,�ۼӱ���
	int	n_trans_total;	//��Һ��������
	int n_bubble_total, n_liquid_total;	//������¼���ݡ�Һ������
	int n_trans[100];	//������¼��Һ����λ������
	{
		n = 1;
		if (sort[n] == 3 || sort[n] == 4)	//Һ������Һ������
		{
			transfer_point_number = 1;
			n_trans[1] = 1;	//��[1]����Һ�����λ��Ϊ1
		}
		else
		{
			transfer_point_number = 0;
		}              /*��ʾ���ݺ�Һ�������������*/
		n = 2;
		while (n <= length)
		{
			if (sort[n] == 3 || sort[n] == 4)
			{
				transfer_point_number = transfer_point_number + 1;
				n_trans[transfer_point_number] = n;                  /*ʹ��n_trans��������¼ÿ���������λ��*/
			}
			n = n + 1;
		}
		n_trans_total = transfer_point_number;                 /*��¼����������*/
		n_bubble_total = n_trans_total / 2;                      /*���ݵ�����*/
		n_liquid_total = n_trans_total / 2;                      /*Һ��������*/
	}
	
    /* ��ʼ������������part01��������Һ���ĳ��ȣ�*/
	n=1;
	while (n<=length)
	{
		if (sort[n]==3)	//Һ-�������壬���һ����Һ���ұ�һ������ 
		{
			x_lleft[n][1]=dx/2;	//�����±�n��ʾ��������ţ�1��ʾʱ�䲽��
			x_lright[n][1]=0;
			x_v[n][1]=dx/2;
		}
		if (sort[n]==4)	//��-Һ�����壬���һ���������ұ�һ����Һ
		{
			x_lleft[n][1]=0;
			x_lright[n][1]=dx/2;
			x_v[n][1]=dx/2;
		}
		if (sort[n]==1)	//Һ������
		{
			x_lleft[n][1]=dx;	//������Ϊɶ����dx/2
			x_lright[n][1]=dx;	//???Ϊɶ����dx/2
			x_v[n][1]=0;
		}
		if (sort[n]==2)	//��������
		{
			x_lleft[n][1]=0;
			x_lright[n][1]=0;
			x_v[n][1]=dx;
		}
		n=n+1;
	}
	/* ��ʼ������������part02��Һ���ٶȣ��¶ȣ�ճ�ȣ����ȣ��ȵ��ʣ�����Ǳ�ȣ�
	������ܶȣ��¶ȣ�ѹ��������Ǳ��,�����¶ȣ�Һ��ѹ����*/
	for (n=1;n<=length;n=n+1)
	{
		if (sort[n]==1)	//Һ������
		{
			v_lleft[n][1] = 0;	//��ʼ�ٶ�Ϊ0
			v_lright[n][1] = 0;
			T_lleft[n][1] = T_c;	//��ʼ�¶����������¶���ͬ
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

		if (sort[n]==2)	//��������
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
									
			T_v[n][1] = T_c;	//��ʼ�¶�
			T_sat[n][1] = T_c;
			ro_v[n][1]= ro_sat_given_T(T_v[n][1]);	//�����ʼ�ܶ�Ϊ����������ܶ�
			P_v[n][1]=p_sat_given_T(T_v[n][1]);	//��ʼѹ��Ϊ����ѹ��
			h_fg_v[n]= h_fg_given_T(T_v[n][1]);	//���ʽ��Һ���һ��
		}

		if (sort[n]==3)	//��Һ����
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
			
			T_v[n][1] = T_c;	//��ʼ�¶�
			T_sat[n][1] = T_c;
			ro_v[n][1] = ro_sat_given_T(T_v[n][1]);	//�����ʼ�ܶ�Ϊ����������ܶ�
			P_v[n][1] = p_sat_given_T(T_v[n][1]);	//��ʼѹ��Ϊ����ѹ��
			h_fg_v[n] = h_fg_given_T(T_v[n][1]);	//���ʽ��Һ���һ��
		}

		if (sort[n]==4)	//������Һ
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

			T_v[n][1] = T_c;	//��ʼ�¶�
			T_sat[n][1] = T_c;
			ro_v[n][1] = ro_sat_given_T(T_v[n][1]);	//�����ʼ�ܶ�Ϊ����������ܶ�
			P_v[n][1] = p_sat_given_T(T_v[n][1]);	//��ʼѹ��Ϊ����ѹ��
			h_fg_v[n] = h_fg_given_T(T_v[n][1]);	//���ʽ��Һ���һ��
		}

		T_w[n][1]=T_c;	//�����¶�Ϊ�������¶�
		P_l[n][2]=0;	//Һ���ʼѹ��Ϊ0Pa
	}
	/* �����γ�ʼ�¶�*/
	for (n=1;n<=length;n=n+1)
	{
		if (heat_sig[n]==1)
		{T_w[n][1]=19.67;}	//�����α����¶ȣ������������ֱ����¶�
	}
	/*��ʼ��Һ�嵥���������ϵ��*/
	double Re_lleft[length+1],Re_lright[length+1],f,Pr;
	for (n=1;n<=length;n=n+1)
	{
		if (sort[n]==1||sort[n]==3||sort[n]==8)	//1��ʾҺ��2��ʾҺ-����8��ʾҺ-��-Һ���������Һ�廻��ϵ��
		{
			Re_lleft[n]=d*ro_l*fabs(v_lleft[n][1])/mu_left[n];
			Pr = Pr_l_given_T(T_lleft[n][1]);
			if (Re_lleft[n]<=2300)
			{
				h_lleft[n]=4.364*lamt_left[n]/d;
			}
			else if (Re_lleft[n]<=10000)
			{
				f=pow(1.82*log10(Re_lleft[n])-1.64,-2);	//Ħ��ϵ��
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
		if (sort[n]==1||sort[n]==4||sort[n]==8)	//1��ʾҺ��4��ʾ��-Һ��8��ʾҺ-��-Һ�������Ҳ�Һ�嵥���������ϵ��
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

	/*******ȷ�������ݵ�״̬*********/
	int n_bubble;	//������ţ��ۼӱ���
	double P_sat;	//���ݵı���ѹ����Pa
	double Tr_v[100];	//���������������kg/s������0��ʾҺĤ������������������
	double V_bubble[100][3];	//���������m3
	double m_bubble[100][3];	//����������kg
	double T_vb[100][3];	//�����¶ȣ�K
	double T_sat_vb[100][3];	//���ݱ����¶ȣ�K����������ѹ�����
	double P_vb[100][3];	//����ѹ��,Pa

	if (sort[1]==2||sort[1]==4)	//��1������������������Һ����ζ�����ݿ�������ԭ�㣬������Ҫ�ҵ������ݵ����
	{
		n_bubble=1;	//��ǰ�������Ϊ1
		while (n_bubble<=n_bubble_total)	//�������С�ڵ�����������
		{
			if (n_bubble==1)	//��1������
			{
				Tr_v[n_bubble]=0;
				n=n_trans[2*n_bubble_total];	//��λ�����һ������λ�ã�����1�����ݵ����
				while (n<=length)	//����1������������ԭ����߸����������������ۼ�
				{
					P_sat=p_sat_given_T(T_w[n][1]);	//�����¶ȶ�Ӧ�ı���ѹ��
					T_sat[n][1]=T_sat_given_p(P_v[n][1]);
					if (P_sat>=P_v[n][1])
						{Tr_v_sig[n]=1;}	//����
					else
						{Tr_v_sig[n]=-1;}	//����
					Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][1]*fabs(T_w[n][1]-T_sat[n][1])/h_fg_v[n];
					n=n+1;
				}
				n=1;	//��λ����1�������壬����1������������ԭ���ұ߸����������������ۼ�
				while (n<=n_trans[2*n_bubble-1])	//n��û����1�����ݵ��յ�
				{
					P_sat=p_sat_given_T(T_w[n][1]);
					T_sat[n][1]= T_sat_given_p(P_v[n][1]);
					if (P_sat>=P_v[n][1])
						{Tr_v_sig[n]=1;}	//����
					else
						{Tr_v_sig[n]=-1;}	//����
					Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][1]*fabs(T_w[n][1]-T_sat[n][1])/h_fg_v[n];
					n=n+1;
				}
				/*��1�����ݵ����=          �յ���������+               �����������+                    +�м������������*/
				V_bubble[n_bubble][1]=A_v*(x_v[n_trans[2*n_bubble-1]][1]+x_v[n_trans[2*n_bubble_total]][1]+(length-1-n_trans[2*n_bubble_total]+n_trans[2*n_bubble-1])*dx);
				m_bubble[n_bubble][1]=ro_v[n_trans[2*n_bubble-1]][1]*V_bubble[n_bubble][1];		//�յ��������ܶ�*���		
			}
			else//n_bubble!=1
			{
				Tr_v[n_bubble]=0;
				n=n_trans[2*n_bubble-2];	//��λ���������
				while (n<=n_trans[2*n_bubble-1])	//������㵽�յ�
				{
					P_sat=p_sat_given_T(T_w[n][1]);
					T_sat[n][1]= T_sat_given_p(P_v[n][1]);
					if (P_sat>=P_v[n][1])
						{Tr_v_sig[n]=1;}	//����
					else
						{Tr_v_sig[n]=-1;}	//����
					Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][1]*fabs(T_w[n][1]-T_sat[n][1])/h_fg_v[n];
					n=n+1;
				}
				V_bubble[n_bubble][1]=A_v*(x_v[n_trans[2*n_bubble-2]][1]+x_v[n_trans[2*n_bubble-1]][1]+(n_trans[2*n_bubble-1]-n_trans[2*n_bubble-2]-1)*dx);
				m_bubble[n_bubble][1]=ro_v[n_trans[2*n_bubble-1]][1]*V_bubble[n_bubble][1];
			}
			T_vb[n_bubble][1]=T_v[n_trans[2*n_bubble-1]][1];	//��¼��n_bubble�����ݵ��¶ȣ����յ�������¶ȱ�ʾ
			T_sat_vb[n_bubble][1] = T_sat[n_trans[2 * n_bubble - 1]][1];	//��¼��n_bubble�����ݵı����¶ȣ����յ�������¶ȱ�ʾ
			P_vb[n_bubble][1]=P_v[n_trans[2*n_bubble-1]][1];	//��¼��n_bubble�����ݵ�ѹ�������յ������ѹ����ʾ
			n_bubble=n_bubble+1;	//���μ����2��3����n_bubble_total�����ݵ����������������������¶ȣ�ѹ��
		}
	}
	else//��1�������岻�������忪ʼ��������Һ�忪ʼ��û�����ݿ�����ԭ��
	{
		n_bubble=1;
		while (n_bubble<=n_bubble_total)
		{
			Tr_v[n_bubble]=0;
			n=n_trans[2*n_bubble-1];	//��λ�������ݵ�*���*λ��
			while (n<=n_trans[2*n_bubble])	//����㵽�յ�
				{
					P_sat=p_sat_given_T(T_w[n][1]);
					T_sat[n][1]= T_sat_given_p(P_v[n][1]);
					if (P_sat>=P_v[n][1])
						{Tr_v_sig[n]=1;}	//����
					else
						{Tr_v_sig[n]=-1;}	//����
					Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][1]*fabs(T_w[n][1]-T_sat[n][1])/h_fg_v[n];
					n=n+1;
				}
			V_bubble[n_bubble][1]=A_v*(x_v[n_trans[2*n_bubble-1]][1]+x_v[n_trans[2*n_bubble]][1]+(n_trans[2*n_bubble]-n_trans[2*n_bubble-1]-1)*dx);
			m_bubble[n_bubble][1]=ro_v[n_trans[2*n_bubble-1]][1]*V_bubble[n_bubble][1];	//������������ܶ�
			T_vb[n_bubble][1]=T_v[n_trans[2*n_bubble-1]][1];	//������������¶ȱ�ʾ�����ݵ��¶�
			T_sat_vb[n_bubble][1] = T_sat[n_trans[2 * n_bubble - 1]][1];	//������������¶ȱ�ʾ�����ݵı����¶�
			P_vb[n_bubble][1]=P_v[n_trans[2*n_bubble-1]][1];	//�����������ѹ����ʾ��ѹ�����¶�
			n_bubble=n_bubble+1;	//���μ����2��3��n_bubble_total�����ݵ����������������������¶ȣ�ѹ��
		}	
	}

	/******Һ���ķ�������******/
	double Tr_vl_lleft[100], Tr_vl_lright[100];	//Һ����(��)������ڴ������ʣ�kg/s������0��ʾҺ��������С
	double Tr_vl[100];	//Һ���ܷ��ڴ������ʣ�kg/s������0��ʾҺ��������С
	double dl_vl_left[100],dl_vl_right[100];	//Һ����(��)�������ڷ����ƶ��ľ��룬m������0��ʾҺ�����ȼ�С
	int Tr_vl_lleft_sig[length+1],Tr_vl_lright_sig[length+1];	//Һ������ǣ�1��ʾ���ڣ�0��ʾû�����
	double Tr_vl_left_each[length+1],Tr_vl_right_each[length+1];	//����������ң�������ڴ�������,kg/s
	int n_liquid;	//Һ�����
	double v_minset=1e-3;	//�ٶ�����10-3 m/s
	/*��ÿ��Һ�������99������ʼ��*/
	for (n=1;n<=99;n=n+1)
	{
		Tr_vl_lleft[n]=0;
		Tr_vl_lright[n]=0;
		Tr_vl[n]=0;
		dl_vl_left[n]=0;
		dl_vl_right[n]=0;
	}
	/*ȷ����������ķ��ڱ��*/
	for (n=1;n<=length;n=n+1)
	{
		Tr_vl_left_each[n]=0;	//ÿ�����������Һ������������ʣ�kg/s
		Tr_vl_right_each[n]=0;	//ÿ���������Ҳ�Һ�������������,kg/s
		Tr_vl_lright_sig[n]=0;	//ÿ�����������Һ������ǳ�ʼ��
		Tr_vl_lleft_sig[n]=0;	//ÿ���������Ҳ�Һ������ǳ�ʼ��
		h_b[n]=0;	//���ڻ���ϵ����ʼ��
		if (heat_sig[n]==1)	//��ʾ������
		{
			if (T_lright[n][1]!=0&&(T_w[n][1]>=T_lright[n][1]))	//�����¶ȸ���Һ���¶�
			{
				Tr_vl_lright_sig[n]=1;	//���ڱ��
			}
			if (T_lleft[n][1]!=0&&(T_w[n][1]>=T_lleft[n][1]))
			{
				Tr_vl_lleft_sig[n]=1;	//���ڱ��
			}
		}
	}
	
	/****����Һ�����ҽ����λ��������ͷ����������������������ƶ�����****/
	/**��1��������Ϊ��������-Һ**/
	if (sort[1]==2||sort[1]==4)	
	{
		n_liquid=1;	//��1��Һ����ʼѭ��
		while (n_liquid<=n_liquid_total)
		{
			/*Һ���������������1*/
			if (n_trans[2*n_liquid-1]>=n_heat_start1&&n_trans[2*n_liquid-1]<=n_heat_end1)
			{
				/*Һ���ҽ��治��������1*/
				if (n_trans[2*n_liquid]>n_heat_end1)
				{
					n=n_trans[2*n_liquid-1];	//��λ��Һ�������
					/*������浽������1���յ㣬�����������ȫ���ۼƵ������*/
					while (n<=n_heat_end1)	
					{
						if (fabs(v_lright[n][1])<=v_minset)	//�����壨�ұ�Һ�壩���ٶ�С���ٶ����ޣ����ڻ���ϵ�����ٶȰ����޼���
						{
							h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);
						}
						else//�������ٶȴ��ڵ������ޣ�����ʵ��ֵ������ڻ���ϵ��
						{
							h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);
						}
						if (T_w[n][1]>=T_lright[n][1])//������ڷ���������
						{
							Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//Һ����������1�ķ����������ȫ���ۼƵ������
							Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//ÿ��������ķ��ڴ�����
							if (sort[n]==1)//Һ������ķ��ڴ���������=��
							{
								Tr_vl_left_each[n]=Tr_vl_right_each[n];	
							}
						}
						n=n+1;
					}
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//Һ����������������
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//Һ���������ƶ�����
					/*Һ���ҽ�����������2����������2����㵽Һ���ҽ�����ڴ������ۼӵ��ҽ���*/
					if ((n_trans[2*n_liquid]>=n_heat_start2)&&(n_trans[2*n_liquid]<=n_heat_end2))
					{
						n=n_heat_start2;	//��λ��������2�����
						while (n<=n_trans[2*n_liquid]-1)
						{
							if (fabs(v_lright[n][1])<=v_minset)//���ٶȼ�����ڻ���ϵ��
							{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lright[n])),0.86);}
							else
							{h_b[n]=3000*h_lright[n]*pow((heat_flux_set/(ro_l*fabs(v_lright[n][1])*h_fg_lright[n])),0.86);}
							if (T_w[n][1]>=T_lright[n][1])//�����������
							{
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ������ۼӵ��ҽ���
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//ÿ���������ҽ�����ڴ�����
								if (sort[n]==1)
							    {Tr_vl_left_each[n]=Tr_vl_right_each[n];}	//Һ��������ڴ���������=��
							}
							n=n+1;
						}
						n=n_trans[2*n_liquid];	//��λ��Һ���ҽ������ڿ�����
						if (fabs(v_lleft[n][1])<=v_minset)
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*v_minset*h_fg_lleft[n])),0.86);}
						else
						{h_b[n]=3000*h_lleft[n]*pow((heat_flux_set/(ro_l*fabs(v_lleft[n][1])*h_fg_lleft[n])),0.86);}
						if (T_w[n][1]>=T_lleft[n][1])
						{
							Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//��ǰ�治֮ͬ������x_lleft[n][1]
							Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//��ǰ�治֮ͬ������x_lleft[n][1]
							if (sort[n]==1)
							{Tr_vl_right_each[n]=Tr_vl_left_each[n];}
						}
						if (n_liquid==n_liquid_total)//Һ���ҽ���ķ����������ۼӵ���һ��Һ��������
						{Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];}
						else
						{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];}
					
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ������ڷ��ڲ������ƶ�����
					}
					/*Һ���ҽ��泬����������2��������2�ķ��������ۼӵ��ҽ���*/
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
						if (n_liquid==n_liquid_total)//Һ���ҽ���ķ����������ۼӵ���һ��Һ��������
						{Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];}
						else
						{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];}
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);
					}
				}
				/*Һ���ҽ�����������1��������浽�ҽ�������п����壬���������ۼӵ�����棬Ȼ��ƽ���ָ���������Һ��*/
				else
				{
					n=n_trans[2*n_liquid-1];	//��λ��Һ�������
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
					n=n_trans[2*n_liquid];	//��λ���ҽ������ڿ�����
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
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//�ָ�n_liquidһ��
					if (n_liquid==n_liquid_total)	//�ָ�n_liquid+1һ��
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lleft[n_liquid]/2;}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lleft[n_liquid]/2;}
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//ֻ�������λ�ƣ�����
				}	
			}
			/*Һ���������������2*/
			else if (n_trans[2*n_liquid-1]>=n_heat_start2&&n_trans[2*n_liquid-1]<=n_heat_end2)
			{
				/*Һ���ҽ��泬����������2*/
				if (n_trans[2*n_liquid]>n_heat_end2)
				{
					n=n_trans[2*n_liquid-1];	//��λ�������
					/*������浽������2���յ㣬�����������ȫ���ۼƵ������*/
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
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//Һ������������
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//Һ���������ƶ�����
				}
				/*Һ���ҽ���û����������2���������浽�ҽ�������п����壬���������ۼӵ�����棬Ȼ��ƽ���ָ���������Һ��������*/
				else
				{
					n=n_trans[2*n_liquid-1];	//��λ�������
					while (n<=(n_trans[2*n_liquid]-1))	//������浽�ҽ�����ߵĿ����壬���������ۼ�
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
					n=n_trans[2*n_liquid];	//��λ���ҽ���
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
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//n_liquid����һ��
					if (n_liquid==n_liquid_total)	//n_liquid+1����һ��
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lleft[n_liquid]/2;}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lleft[n_liquid]/2;}
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//λ��ȫ������������
				}
			}
			/*�������1��������1�����֮��*/
			else if (n_trans[2*n_liquid-1]>=1&&n_trans[2*n_liquid-1]<n_heat_start1)
			{
				/*�ҽ�����������1*/
				if (n_trans[2*n_liquid]>=n_heat_start1&&n_trans[2*n_liquid]<=n_heat_end1)
				{
					n=n_heat_start1;	//��λ��������1���
					/*��������1��㵽Һ���ҽ��棬���������ۼӵ��ҽ���*/
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
					n=n_trans[2*n_liquid];	//��λ���ҽ���
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
					/*�ҽ���ķ��������������һ��Һ������*/
					if (n_liquid==n_liquid_total)
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];}
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);
				}
				/*�ҽ�����������1֮����������2֮ǰ*/
				else if (n_trans[2*n_liquid]>n_heat_end1&&n_trans[2*n_liquid]<n_heat_start2)
				{
					n=n_heat_start1;	//��λ��������1���
					/*������1�ķ��������ۼӵ�����棬Ȼ��ƽ���ָ���������Һ��*/
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
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//n_liquid����һ��
					if (n_liquid==n_liquid_total)	//n_liquid+1����һ��
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lleft[n_liquid]/2;}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lleft[n_liquid]/2;}
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//�������ƶ�����
				}
				/*�ҽ�����������2*/
				else if (n_trans[2*n_liquid]>=n_heat_start2&&n_trans[2*n_liquid]<=n_heat_end2)
				{
					n=n_heat_start1;	//��λ��������1���
					/*������1�ķ��������ۼӵ������*/
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
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//�������������n_liquid
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
					n=n_heat_start2;	//��λ��������2���
					/*������2��㵽�ҽ���ķ��������ۼӵ��ҽ���*/
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
					if (n_liquid==n_liquid_total)	//���������������һ��Һ��
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];}
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
				}
				/*�ҽ��泬��������2*/
				else if (n_trans[2*n_liquid]>n_heat_end2)
				{
					n=n_heat_start1;	//��λ��������1���
					/*������1���������ۼӵ������*/
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
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//�������������n_liquid
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
					n=n_heat_start2;	//��λ��������2���
					/*������2���������ۼӵ��ҽ���*/
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
					if (n_liquid==n_liquid_total)	//�������������n_liquid+1
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];}
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
				}

			}
			/*�������������1֮����������2֮ǰ*/
			else if (n_trans[2*n_liquid-1]>n_heat_end1&&n_trans[2*n_liquid-1]<n_heat_start2)
			{
				/*�ҽ�����������2*/
				if (n_trans[2*n_liquid]>=n_heat_start2&&n_trans[2*n_liquid]<=n_heat_end2)
				{
					n=n_heat_start2;	//��λ��������2���
					/*������2��㵽�ҽ���ķ��������ۼӵ��ҽ���*/
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
					if (n_liquid==n_liquid_total)	//�����������䵽��һ��Һ��
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];}
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
				}
				/*�ҽ��泬��������2*/
				else if (n_trans[2*n_liquid]>n_heat_end2)
				{
					n=n_heat_start2;	//��λ��������2���
					/*������2�ķ��������ۼӵ��ҽ��棬ƽ���������������Һ��*/
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
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid]/2;	//n_liquid����һ��
					if (n_liquid==n_liquid_total)	//n_liquid+1����һ��
					{Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid]/2;}
					else
					{Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid]/2;}
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ�����ƶ�����
				}
			}
			/*Һ�����ҽ���λ�÷���������ϣ�������һ��Һ��*/
			n_liquid=n_liquid+1;
		}
	}
	/**��1��������ΪҺsort[1]==1����Һ-��sort[1]==3����Һ-��-Һsort[1]==8**/
	else                                
	{
		n_liquid=1;   //��1��Һ���������ǣ�����Ϊ����ԭ�㽫��1��Һ���ֳ������Σ��ȿ���ԭ���ұ���һ�Σ��ٿ���ԭ�������һ��
		/*��1��Һ�����ҽ��棨n_trans[1]����������1*/
		if (n_trans[1]>=n_heat_start1&&n_trans[1]<=n_heat_end1)
		{
			n=n_heat_start1;	//��λ��������1���
			/*������1��㵽�ҽ���ķ��������ۼӵ��ҽ���*/
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
			Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�������������n_liquid
			dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
		}
		/*��1��Һ���ҽ�����������1֮��������2֮ǰ*/
		else if (n_trans[1]>n_heat_end1&&n_trans[1]<n_heat_start2)
		{
			n=n_heat_start1;	//��λ��������1���
			/*������1���������ۼӵ��ҽ��棬ƽ�������������Һ��*/
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
			Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid]/2;	//n_liquid����һ��
			Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lright[n_liquid]/2;	//n_liquid_total����һ��
			dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
		}
		/*��1��Һ���ҽ�����������2*/
		else if (n_trans[1]>=n_heat_start2&&n_trans[1]<=n_heat_end2)
		{
			n=n_heat_start1;	//��λ��������1���
			/*������1�ķ��������ۼӵ�����棿����*/
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
			Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid];	//�������������n_liquid_total
			dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
			n=n_heat_start2;	//��λ��������2���
			/*������2��㵽�ҽ���ķ��������ۼӵ��ҽ���*/
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
			Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�������������n_liquid
			dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
		}

		/*��1��Һ�������(n_trans[2*n_liquid_total]����������1*/
		if (n_trans[2*n_liquid_total]>=n_heat_start1&&n_trans[2*n_liquid_total]<=n_heat_end1)
		{
			n=n_trans[2*n_liquid_total];	//��λ�������
			/*����浽������1�յ���������ۼӵ������*/
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
			Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid];	//�������������n_liquid_total
			dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
			n=n_heat_start2;	//��λ��������2���
			/*������2���������ۼӵ��ҽ���*/
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
			Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�������������n_liquid=1
			dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
		}
		/*��1��Һ���������������1֮����������2֮ǰ*/
		else if (n_trans[2*n_liquid_total]>n_heat_end1&&n_trans[2*n_liquid_total]<n_heat_start2)
		{
			n=n_heat_start2;	//��λ��������2���
			/*������2�ķ��������ۼӵ�����棬ƽ�������������Һ��*/
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
			Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid]/2;	//�����n_liquid_totalһ��
			Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//�����n_liquidһ��
			dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
		}
		/*��1��Һ���������������2*/
		else if (n_trans[2*n_liquid_total]>=n_heat_start2&&n_trans[2*n_liquid_total]<=n_heat_end2)
		{
			n=n_trans[2*n_liquid_total];	//��λ�������
			/*����浽������2�յ�ķ��������ۼӵ������*/
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
			Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid];	//�������������n_liquid_total
			dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
		}

		n_liquid=2;	//��2��Һ����ʼѭ��
		while (n_liquid<=n_liquid_total)
		{
			/*�����n_trans[2*n_liquid-2]��������1*/
			if (n_trans[2*n_liquid-2]>=n_heat_start1&&n_trans[2*n_liquid-2]<=n_heat_end1)
			{
				/*�ҽ��泬��������1*/
				if (n_trans[2*n_liquid-1]>n_heat_end1)
				{
					n=n_trans[2*n_liquid-2];	//��λ�������
					/*����浽������1�յ�ķ��������ۼӵ������*/
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
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//�������������n_liquid-1
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
					/*�ҽ�����������2*/
					if ((n_trans[2*n_liquid-1]>=n_heat_start2)&&(n_trans[2*n_liquid-1]<=n_heat_end2))
					{
						n=n_heat_start2;	//��λ��������2���
						/*������2��㵽�ҽ���ķ��������ۼӵ��ҽ���*/
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
						n=n_trans[2*n_liquid-1];	//��λ���ҽ�������壬������������������ۼӵ��ҽ���
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
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�������������n_liquid
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
					}
					/*�ҽ��泬��������2*/
					if (n_trans[2*n_liquid-1]>n_heat_end2)
					{
						n=n_heat_start2;	//��λ��������2���
						/*������2���������ۼӵ��ҽ���*/
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
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�������������n_liquid
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
					}
				}
				/*�ҽ�����������1*/
				else
				{
					n=n_trans[2*n_liquid-2];	//��λ�������
					/*������浽�ҽ�����������ۼӵ�����棬ƽ�����䵽������������*/
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
					n=n_trans[2*n_liquid-1];	//��λ���ҽ��棬������������������ۼӵ������
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
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid]/2;	//����һ���n_liquid-1
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//����һ���n_liquid
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
				}	
			}
			/*�������������2*/
			else if (n_trans[2*n_liquid-2]>=n_heat_start2&&n_trans[2*n_liquid-2]<=n_heat_end2)
			{
				/*�ҽ��泬��������2*/
				if (n_trans[2*n_liquid-1]>n_heat_end2)
				{
					n=n_trans[2*n_liquid-2];	//��λ�������
					/*����浽������2�յ�ķ��������ۼӵ������*/
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
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//�������������n_liquid-1
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
				}
				/*�ҽ�����������2*/
				else
				{
					n=n_trans[2*n_liquid-2];	//��λ�������
					/*������浽�ҽ�����������ۼӵ�����棬ƽ�����䵽������������*/
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
					n=n_trans[2*n_liquid-1];	//��λ���ҽ�������壬������������������ۼӵ������
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
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid]/2;	//����һ���n_liquid-1
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//����һ���n_liquid
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
				}
			}
			/*�������1��������1���֮��*/
			else if (n_trans[2*n_liquid-2]>=1&&n_trans[2*n_liquid-2]<n_heat_start1)
			{
				/*�ҽ�����������1*/
				if (n_trans[2*n_liquid-1]>=n_heat_start1&&n_trans[2*n_liquid-1]<=n_heat_end1)
				{
					n=n_heat_start1;	//��λ��������1���
					/*������1��㵽�ҽ�����������ۼӵ��ҽ���*/
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
					n=n_trans[2*n_liquid-1];	//��λ���ҽ�������壬��������������������ۼӵ��ҽ���
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
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�����n_liquid
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
				}
				/*�ҽ�����������1֮��������2֮ǰ*/
				else if (n_trans[2*n_liquid-1]>n_heat_end1&&n_trans[2*n_liquid-1]<n_heat_start2)
				{
					n=n_heat_start1;	//��λ��������1���
					/*������1�ķ��������ۼӵ�����棬ƽ�������������������*/
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
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid]/2;	//����һ���n_liquid-1�����������
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//����һ���n_liquid�����Ҳ�����
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//�������ƶ�����
				}
				/*�ҽ�����������2*/
				else if (n_trans[2*n_liquid-1]>=n_heat_start2&&n_trans[2*n_liquid-1]<=n_heat_end2)
				{
					n=n_heat_start1;	//��λ��������1���
					/*������1���������ۼӵ������*/
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
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//�����n_liquid-1�����������
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
					n=n_heat_start2;	//��λ��������2���
					/*������2��㵽�ҽ�����������ۼӵ��ҽ���*/
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
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�������������n_liquid�����Ҳ�����
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
				}
				/*�ҽ��泬��������2*/
				else if (n_trans[2*n_liquid-1]>n_heat_end2)
				{
					n=n_heat_start1;	//��λ��������1���
					/*������1���������ۼӵ������*/
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
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//�������������n_liquid-1�����������
					dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
					n=n_heat_start2;	//��λ��������2���
					/*������2���������ۼӵ��ҽ���*/
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
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�������������n_liquid�����Ҳ�����
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
				}
			}
			/*�������������1֮��������2֮ǰ*/
			else if (n_trans[2*n_liquid-2]>n_heat_end1&&n_trans[2*n_liquid-2]<n_heat_start2)
			{
				/*�ҽ�����������2*/
				if (n_trans[2*n_liquid-1]>=n_heat_start2&&n_trans[2*n_liquid-1]<=n_heat_end2)
				{
					n=n_heat_start2;	//��λ��������2���
					/*������2���������ۼӵ��ҽ���*/
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
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�������������n_liquid�����Ҳ�����
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
				}
				/*�ҽ��泬��������2*/
				else if (n_trans[2*n_liquid-1]>n_heat_end2)
				{
					n=n_heat_start2;	//��λ��������2���
					/*������2���������ۼӵ��ҽ��棬ƽ�������������������*/
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
					Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lright[n_liquid]/2;	//�����ұ߸�n_liquid-1�����������
					Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid]/2;	//����һ���n_liquid�����Ҳ�����
					dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
				}
			}
			/*����������ϣ���һ��Һ��*/
			n_liquid=n_liquid+1;
		}

	}

	printf("��ʼ�����\n");	//��Ļ���

	FILE *fid0;	//����һ����fid0���ļ�ָ��
	fid0 =fopen("state.txt","a+");	//����a+ģʽ��һ������state.txt���ļ�������state.txt���ļ�ָ���fid0,�쳣�򷵻�NULL
	/*a+��ʾ�Ը��ӷ�ʽ�򿪿ɶ�д���ļ������ļ������ڣ���Ὠ�����ļ�������ļ����ڣ�д������ݻᱻ�ӵ��ļ�β��
	���ļ�ԭ�ȵ����ݻᱻ������ԭ����EOF��������*/
	fprintf(fid0,"��ʼ�����\n");	//������ʼ�����\n��д����fid0ָ�����ļ�����ʽ��printf���ƣ�����ָ��д����������
	fclose(fid0);	//�ر���fid0ָ�����ļ�,���ز��������0��EOF

	FILE *fid1=fopen("ro_v.txt","a+");	//����a+ģʽ��һ������ro_v.txt���ļ�������ro_v.txt���ļ�ָ���fid1
	FILE *fid2=fopen("T_w.txt","a+");	//����a+ģʽ��һ������T_w.txt���ļ�������T_w.txt���ļ�ָ���fid2
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
	
	/*д���һ�б�ͷ��ʾn��1��length��Ȼ����*/
	int location = 1;       //���ʱ���������꣬�൱��n*/
	for (location = 1; location <= length; location++)
	{
		fprintf(fid1, " %d", location);	//�ո����ڷ��У�n=1�ڵ�2��
		fprintf(fid2, " %d", location);	//�ո����ڷ��У�n=1�ڵ�2��
		fprintf(fid4, " %d", location);	//�ո����ڷ��У�n=1�ڵ�2��
		fprintf(fid6, " %d", location);	//�ո����ڷ��У�n=1�ڵ�2��
		fprintf(fid8, " %d", location);	//�ո����ڷ��У�n=1�ڵ�2��
		fprintf(fid9, " %d", location);	//�ո����ڷ��У�n=1�ڵ�2��
		fprintf(fid10, " %d", location);	//�ո����ڷ��У�n=1�ڵ�2��
		fprintf(fid16, " %d", location);	//�ո����ڷ��У�n=1�ڵ�2��
		fprintf(fid20, " %d", location);	//�ո����ڷ��У�n=1�ڵ�2��
		fprintf(fid21, " %d", location);	//�ո����ڷ��У�n=1�ڵ�2��
		fprintf(fid22, " %d", location);	//�ո����ڷ��У�n=1�ڵ�2��
	}
	/*��ͷ����n֮����*/
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

	fclose(fid1);	//�ر���fid1ָ�����ļ�
	fclose(fid2);	//�ر���fid2ָ�����ļ�
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

	/*��ʼ����һʱ��ڵ��ֵ*/
	for (n=1;n<=length;n=n+1)
	{
		x_lleft[n][2]=0;	//���������Һ������
		x_lright[n][2]=0;	//�������ұ�Һ������
		v_lleft[n][2]=0;	//���������Һ���ٶ�
		v_lright[n][2]=0;
		T_lleft[n][2]=0;	//���������Һ���¶�
		T_lright[n][2]=0;
		x_v[n][2]=0;	//����������������
		ro_v[n][2]=0;	//�����ܶ�
		T_v[n][2]=0;	//�����¶�
		T_sat[n][2] = 0;	//���������¶�
		P_v[n][2]=0;	//����ѹ��
	}

	/*******��ʼ��ѭ��*********/

    int i=1;	//i��ʱ��ڵ�ѭ�����Ʊ���
	int kk=0;	//�������������
	double sensible_c_total=0;	//�����������ȴ����������ʱ�䲽���µ��ۼ�������W
	double sensible_e_total=0;	//�����������ȴ����������ʱ�䲽���µ��ۼ�������W
	double latent_c_total=0;	//��������Ǳ�ȴ����������ʱ�䲽���µ��ۼ�������W
	double latent_e_total=0;	//��������Ǳ�ȴ����������ʱ�䲽���µ��ۼ�������W

	while (i<=nn)
	{
		int cross_positive=0;	//Һ��������ƶ���ǣ�=1��ʾ������ԭ�㣬��1�����ݱ��2�� 
		int cross_negative=0;	//Һ��������ƶ���ǣ�=1��ʾ������ԭ�㣬��1�����ݱ����һ��
		double sensible_c=0;	//���������ȣ�Һ�嵥��������ȣ���W
		double latent_c=0;	//������Ǳ��,W
		double sensible_e=0;	//����������,W
		double latent_e=0;	//������Ǳ��,W

		/*******���㴫���������Ⱥ�Ǳ�ȣ�*********/
		n=1;
		while (n<=length)
		{
			if (heat_sig[n]==-1)	//-1��ʾ����
			{
				/*�������������Ⱥ�Ǳ��*/
				if (sort[n]==1)	//1��ʾҺ������
				{
					sensible_c=sensible_c+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d;	//�����������������ߵ�Һ���¶ȳ��Ȼ���ϵ������ʾ
				}
				else //������������Ŀ����壬2����3Һ-����4��-Һ��8Һ-��-Һ
				{
					sensible_c=sensible_c+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d;	//������������Һ���¶ȳ��Ȼ���ϵ��
					sensible_c=sensible_c+h_lright[n]*(T_w[n][1]-T_lright[n][1])*x_lright[n][1]*pi*d;	//����������ұ�Һ���¶ȳ��Ȼ���ϵ��
					latent_c=latent_c+Tr_v_sig[n]*h_v[n]*fabs(T_w[n][1]-T_sat[n][1])*x_v[n][1]*pi*d;	//ҺĤ������/��������
				}
			}
			if (heat_sig[n]==1)	//1��ʾ������
			{
				/*�������������Ⱥ�Ǳ��*/
				if (sort[n]==1)		//Һ������			
				{
					sensible_e=sensible_e+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d;	//������������Һ�����
					latent_e=latent_e+Tr_vl_lleft_sig[n]*h_b[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d;	//���ڣ����Һ�����
				}
				else//����������
				{
					sensible_e=sensible_e+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d;	//�����������
					sensible_e=sensible_e+h_lright[n]*(T_w[n][1]-T_lright[n][1])*x_lright[n][1]*pi*d;	//�����������
					latent_e=latent_e+Tr_v_sig[n]*h_v[n]*fabs(T_w[n][1]-T_sat[n][1])*x_v[n][1]*pi*d
							+Tr_vl_lleft_sig[n]*h_b[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
							+Tr_vl_lright_sig[n]*h_b[n]*(T_w[n][1]-T_lright[n][1])*x_lright[n][1]*pi*d;	//ҺĤ����+�����+�ҷ���
				}
			}
			n=n+1;
		}

		/*******��������¶ȣ���1�������1�������嵥�����ǣ�*******/
		n=1;
		if (plain[n]==0)	//0��ʾ���¿ɱ�
		{
			if (sort[n]==1)	//Һ������
			{
				T_w[n][2]=T_w[n][1]+dt/(c_w*ro_w*A_w*dx)*(heat_flux[n]*dx-h_c[n]*pi*d_out*(T_w[n][1]-T_out)*dx
						-Tr_v_sig[n]*h_v[n]*fabs(T_w[n][1]-T_sat[n][1])*x_v[n][1]*pi*d  //��һ����ࣿ����
						-Tr_vl_lleft_sig[n]*h_b[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
						-h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*x_lleft[n][1]*pi*d
						-lamt_w*A_w*(T_w[n][1]-T_w[length][1])/dx-lamt_w*A_w*(T_w[n][1]-T_w[n+1][1])/dx);	//������������
			}
			else//����������
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
		else//�����ѱ��趨����Ϊ��ʼֵ��ֻ��������
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

		/*******��ʼ��ÿ��������������*******/
		for (n=1;n<=length;n=n+1)
			{Tr_v_sig[n]=0;}//0��ʾ�����

		/*******����һЩ��*******/
		int jump;	//������ƶ���ǣ�0��ʾû���Ƴ�ԭ�����壬1��ʾ�����Ƴ�ԭ������
		double dm_left[100], dm_right[100];	//����������䵼�µ�Һ�����ң����澻��������,����0��ʾҺ����������
		double dl_left[100], dl_right[100];	//����������䵼�µ�Һ������)�������ƶ����룬����0��ʾҺ����������
		double G_force[100];	//Һ��������N/m2
		double l[100];	//Һ�����ȣ�
		double C_ll;	//Һ��������ϵ��������0��ʾ���ң�С��0��ʾ����
		double K;	//��ͷ����ϵ��

		/*******������ʼ�㣨n=1��״̬�ķ������ۣ���һ���֣���������-Һ����ʱҺ������������������*********/
		/*******�����Һ���������仯���������á�Һ���ܳ��ȡ��ٶȡ��¶�*********/
		if (sort[1]==2||sort[1]==4)
		{
			n_liquid=1;
			while (n_liquid<=n_liquid_total)
			{
				jump=0;
				/*Һ������澻������������Ϊ���룬��Ϊ������=������������仯��һ��*/
				dm_left[n_liquid]=-Tr_v[n_liquid]/2*dt;		
				dl_left[n_liquid]=dm_left[n_liquid]/A_l/ro_l;	//Һ��������ƶ�����
				/*Һ���ҽ��澻������������Ϊ���룬��Ϊ������=�Ҳ����������仯��һ��*/
				if (n_liquid==n_liquid_total)
					{dm_right[n_liquid]=-Tr_v[1]/2*dt;}
				else
					{dm_right[n_liquid]=-Tr_v[n_liquid+1]/2*dt;}
				dl_right[n_liquid]=dm_right[n_liquid]/A_l/ro_l;    //Һ���ҽ����ƶ�����
				/*���������������ϵ��������ϵ����Һ�峤��*/
				n=n_trans[2*n_liquid-1];	//��λ��Һ�������
				G_force[n_liquid]=0;	//������ʼ��Ϊ0
				l[n_liquid]=0;
				C_ll=0;
				K=0;
				G_force[n_liquid]=G_force[n_liquid]+ro_l*g*c_gravity[n]*x_lright[n][1];	//��ÿ��������������ۼӣ�x_lright[n][1]
				l[n_liquid]=l[n_liquid]+x_lright[n][1];	//��ÿ���������Һ�峤���ۼӣ�x_lright[n][1]
				Re_lright[n]=d*ro_l*fabs(v_lright[n][1])/mu_right[n];	//��������ŵ����v_lright[n][1]
				if (Re_lright[n]==0)	//ȷ������ϵ����Ĭ��Ϊ��ֵ
					{C_l[n]=0;}
				else if (Re_lright[n]<=1180)
					{C_l[n]=16/Re_lright[n];}
				else
					{C_l[n]=0.078*pow(Re_lright[n],-0.2);}
				if((Re_lright[n]>0)&&(n%210>=101)&&(n%210<=210))	//��ͷ����������ʧ
					{K=1000/Re_lright[n]+0.1*(1+4/pow(2.3/25.4,0.3));
					C_l[n]=C_l[n]+K/4*2.3/110;}
				if (v_lright[n_trans[2*n_liquid-1]][1]>0)	//������ٶȴ���0�������������˶�
					{C_l[n]=-C_l[n];}	//����ϵ��Ϊ����������ָ�򸺷���
				C_ll=C_ll+C_l[n]*x_lright[n][1];	//����������ϵ������Һ�峤���ۼӵ�������ϵ��
				/*��������ҵ����п����壬����ϵ��������ϵ����Һ������*/
				n=n+1;
				while (n<=n_trans[2*n_liquid])
				{
					G_force[n_liquid]=G_force[n_liquid]+ro_l*g*c_gravity[n]*x_lleft[n][1];	//��ÿ��������������ۼӣ�x_lleft[n][1]         
					l[n_liquid]=l[n_liquid]+x_lleft[n][1];	//��ÿ���������Һ�峤���ۼӣ�x_lleft[n][1]
					Re_lleft[n]=d*ro_l*fabs(v_lleft[n][1])/mu_left[n];	//��������ŵ����v_lleft[n][1]
					if (Re_lleft[n]==0)
						{C_l[n]=0;}
					else if (Re_lleft[n]<=1180)
						{C_l[n]=16/Re_lleft[n];}
					else
						{C_l[n]=0.078*pow(Re_lleft[n],-0.2);}
					if((Re_lleft[n]>0)&&(n%210>=101)&&(n%210<=210))	//��ͷ����������ʧ
						{K=1000/Re_lleft[n]+0.1*(1+4/pow(2.3/25.4,0.3));
						C_l[n]=C_l[n]+K/4*2.3/110;}
					if (v_lleft[n_trans[2*n_liquid-1]][1]>0)
						{C_l[n]=-C_l[n];}
					C_ll=C_ll+C_l[n]*x_lleft[n][1];
					n=n+1;
				}
				C_ll=C_ll/l[n_liquid];	//����������Һ�峤�ȼ�Ȩƽ����������ϵ��

				/**Һ���������̣����߶�������Һ�������A_l���õ�������������ٶ�**/
				v_lright[n_trans[2*n_liquid-1]][2]=1/(l[n_liquid]+dl_left[n_liquid]+dl_right[n_liquid]-dl_vl_left[n_liquid]-dl_vl_right[n_liquid])
						*(v_lright[n_trans[2*n_liquid-1]][1]*l[n_liquid]
						+dt/ro_l*(G_force[n_liquid]+P_v[n_trans[2*n_liquid-1]][1]-P_v[n_trans[2*n_liquid]][1]+2*C_ll*ro_l/d*l[n_liquid]*pow(v_lright[n_trans[2*n_liquid-1]][1],2)));

				/*�������������v_lֵ����������ͬһҺ���Ŀ�����*/
				n=n_trans[2*n_liquid-1]+1;	//��������ұߵĿ����忪ʼ
				while (n<=n_trans[2*n_liquid]-1)
				{
					v_lright[n][2]=v_lright[n_trans[2*n_liquid-1]][2];
					v_lleft[n][2]=v_lright[n_trans[2*n_liquid-1]][2]; 
					n=n+1;
				}
				n=n_trans[2*n_liquid];	//�ҽ��浥�����ǣ�������Ҫv_lright[n][2]
				v_lleft[n][2]=v_lright[n_trans[2*n_liquid-1]][2];
		
				/*�����ж��Ƿ���λ�õ�ƫ�ƣ�Һ���������յ������λ�øı䣩*/
				/*Һ���������*/
				n=n_trans[2*n_liquid-1];
				x_lright[n][2]=x_lright[n][1]-v_lright[n][1]*dt+dl_left[n_liquid]-dl_vl_left[n_liquid];	//������µ�Һ�峤��
				if (sort[n]!=8)	//����Һ-��-Һ
					{x_lleft[n][2]=x_lleft[n][1];}	//���Һ�峤�ȱ��ֲ���
				/*��������Χ�����һ����Խһ�����ϵĿ����壩������ѭ��while (n_liquid<=n_liquid_total)*/
				if (x_lright[n][2]<=-dx||x_lright[n][2]>=2*dx)	
				{
					codenumber=1;	//��¼�쳣
					break;	//����ѭ��
				}
				/*�����������ƫ�ƣ����������ԭ��������*/
				if (x_lright[n][2]>0&&x_lright[n][2]<dx)	
				{
					/*�����Һ�����������*/
					T_lright[n][2]=1/x_lright[n][2]*(x_lright[n][1]*T_lright[n][1]+dt/(c_pl_right[n]*ro_l*A_l)
						*(-v_lright[n][1]*A_l*ro_l*c_pl_right[n]*(x_lleft[n+1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[n+1][1])/(x_lright[n][1]+x_lleft[n+1][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]
						-2*A_l*(T_lright[n][1]-T_lleft[n+1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[n+1][1]/lamt_left[n+1])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));	//Һ�崢��=����������+���洫��+������+��������-����
					if (T_lright[n][2]>=T_max)
					{
						T_lright[n][2]=T_max;
					}
					if (T_lright[n][2]<=T_min)
					{
						T_lright[n][2]=T_min;
					}
					x_v[n][2]=dx-x_lright[n][2]-x_lleft[n][2];	//������µ����峤��
				}
				/*����������ƶ�һ�������壬��n��n+1�ϲ�Ϊһ��������*/
				else if (x_lright[n][2]<=0)
				{
					x_lright[n+1][2]=dx+x_lright[n][2];	//������µ�Һ�峤��
                    x_v[n+1][2]=dx-x_lright[n+1][2];	//������µ����峤��
                    x_lleft[n+1][2]=0;
                    T_lleft[n+1][2]=0;
                    v_lleft[n+1][2]=0;
                    sort[n+1]=4;	//4��ʾ��-Һ
                    jump=1;	//��¼���������1��������
                      
                    /*�����Һ�����������*/
					T_lright[n+1][2]=1/(x_lright[n+1][2]*c_pl_right[n+1])*(x_lright[n][1]*T_lright[n][1]*c_pl_right[n]+x_lleft[n+1][1]*T_lleft[n+1][1]*c_pl_left[n+1]
						+dt/(A_l*ro_l)*(-v_lleft[n+1][1]*A_l*ro_l*c_pl_right[n+1]*(x_lleft[n+2][1]*T_lright[n+1][1]+x_lright[n+1][1]*T_lleft[n+2][1])/(x_lright[n+1][1]+x_lleft[n+2][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]+h_lleft[n+1]*(T_w[n+1][1]-T_lleft[n+1][1])*pi*d*x_lleft[n+1][1]
						-2*A_l*(T_lright[n+1][1]-T_lleft[n+2][1])/(x_lright[n+1][1]/lamt_right[n+1]+x_lleft[n+2][1]/lamt_left[n+2])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));	//Һ�崢��=����������+���洫��+������+��������-���ڣ��������������������δ����x_lright[n][1]???
					if (T_lright[n+1][2]>=T_max)
					{
						T_lright[n+1][2]=T_max;
					}
					if (T_lright[n+1][2]<=T_min)
					{
						T_lright[n+1][2]=T_min;
					}

					if (sort[n]!=8)	//����Һ-��-Һ������-Һ
					{
						x_lright[n][2]=0;
                        x_lleft[n][2]=0;
                        x_v[n][2]=dx;
                        T_lleft[n][2]=0;
                        T_lright[n][2]=0;
                        v_lleft[n][2]=0;
                        v_lright[n][2]=0;
                        sort[n]=2;	//����������ƶ��������
					}
					else  //nΪҺ-��-Һ
					{
						x_lright[n][2]=0;
                        x_v[n][2]=dx-x_lleft[n][2];
                        sort[n]=3;	//����������ƶ������Һ-��
                        T_lright[n][2]=0;
                        v_lright[n][2]=0;
					}
				}
				/*x_lright[n][2]>=dx������������ƶ�һ�������壬��n��n-1�ϲ�Ϊһ��������*/
				else 
				{
					/*�����Һ�����������*/
					T_lright[n][2]=1/x_lright[n][2]*(x_lright[n][1]*T_lright[n][1]+dt/(c_pl_right[n]*ro_l*A_l)
						*(-v_lright[n][1]*A_l*ro_l*c_pl_right[n]*(x_lleft[n+1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[n+1][1])/(x_lright[n][1]+x_lleft[n+1][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]
						-2*A_l*(T_lright[n][1]-T_lleft[n+1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[n+1][1]/lamt_left[n+1])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));	//���������ԭ�������ڵķ���һ��
					if (T_lright[n][2]>=T_max)
					{
						T_lright[n][2]=T_max;
					}
					if (T_lright[n][2]<=T_min)
					{
						T_lright[n][2]=T_min;
					}

					if (n==1)	//�����ԭ���ڵ�1�������壬�����ƶ������һ��������
					{
						if (sort[length]==3)	//���һ��������ΪҺ-��
						{
							T_lright[length][2]=T_lright[n][2];
                            v_lright[length][2]=v_lright[n][2];
                            x_lright[length][2]=x_lright[n][2]-dx;
                            x_v[length][2]=dx-x_lright[length][2]-x_lleft[length][2];
                            sort[length]=8;	//����������ƶ����Һ-��-Һ
						}
						else  //���һ��������Ϊ��
						{
							x_lright[length][2]=x_lright[n][2]-dx;
                            x_lleft[length][2]=0;
                            //x_lright[length][2]=x_lright[n][2]-dx;
                            T_lright[length][2]=T_lright[n][2];
                            T_lleft[length][2]=0;
                            v_lright[length][2]=v_lright[n][2];
                            v_lleft[length][2]=0;
                            sort[length]=4;	//����������ƶ������-Һ
						}
						cross_negative=1;	//��ʾҺ������渺����ԭ�㣬������ŷ����仯
					}
					else  //�����ԭ�����ڵ�1��������
					{
						if (sort[n-1]==3)	//Һ-��
						{
							T_lright[n-1][2]=T_lright[n][2];
                            v_lright[n-1][2]=v_lright[n][2];
                            x_lright[n-1][2]=x_lright[n][2]-dx;
                            x_v[n-1][2]=dx-x_lright[n-1][2]-x_lleft[n-1][2];
                            sort[n-1]=8;	//���Һ-��-Һ
						}
						else  //��
						{
							x_lright[n-1][2]=x_lright[n][2]-dx;
                            x_lleft[n-1][2]=0;
                            //x_lright[n-1][2]=x_lright[n][2]-dx;
                            T_lright[n-1][2]=T_lright[n][2];
                            T_lleft[n-1][2]=0;
                            v_lright[n-1][2]=v_lright[n][2];
                            v_lleft[n-1][2]=0;
                            sort[n-1]=4;	//�����-Һ
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
                    sort[n]=1;	//ԭ������������Һ

				}

				/*Һ�����м䲿��*/
				if (jump==1)	//����������ƶ�һ��������
				{
					n=n_trans[2*n_liquid-1]+2;	//��λ����������Ҳ�Ŀ�����
					while (n<=n_trans[2*n_liquid]-1)
					{
						/*Һ���������������*/
						T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
							*((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
							-(x_lleft[n+1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[n+1][1])/(x_lright[n][1]+x_lleft[n+1][1]))
							+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
							-2*A_l*(T_lright[n][1]-T_lleft[n+1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[n+1][1]/lamt_left[n+1]));	//Һ�崢�������=����������+���洫��+������
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
				else  //�������ԭ������������ƶ�һ��������
				{
					n=n_trans[2*n_liquid-1]+1;	//��λ��������Ҳ�Ŀ�����
					while (n<=n_trans[2*n_liquid]-1)
					{
						/*Һ���������������*/
						T_lleft[n][2]=T_lleft[n][1]+dt/(dx*c_pl_left[n]*ro_l*A_l)*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]
							*((x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lright[n-1][1]+x_lleft[n][1])
							-(x_lleft[n+1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[n+1][1])/(x_lright[n][1]+x_lleft[n+1][1]))
							+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
							-2*A_l*(T_lright[n][1]-T_lleft[n+1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[n+1][1]/lamt_left[n+1]));	//ͬ��
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
				jump=0;	//���ã�����Ϊʲô������

				/*Һ�����ҽ���*/
				n=n_trans[2*n_liquid];
				x_lleft[n][2]=x_lleft[n][1]+v_lleft[n][1]*dt+dl_right[n_liquid]-dl_vl_right[n_liquid];	//�ҽ��������Һ�峤��
				if (sort[n]!=8)	//����Һ-��-Һ����ֻ����Һ-��
					{x_lright[n][2]=x_lright[n][1];}
				/*��������Χ�����һ����Խһ�����ϵĿ����壩������ѭ��while (n_liquid<=n_liquid_total)*/
				if (x_lleft[n][2]<=-dx||x_lleft[n][2]>=2*dx)
				{
					codenumber=1;	//��¼�쳣
					break;	//����ѭ��
				}
				/*�����������ƫ�ƣ��ҽ�������ԭ��������*/
				if (x_lleft[n][2]>0&&x_lleft[n][2]<dx)
				{
					/*Һ���ҽ������������*/
					T_lleft[n][2]=1/x_lleft[n][2]*(x_lleft[n][1]*T_lleft[n][1]+dt/(c_pl_left[n]*ro_l*A_l)
						*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]*(x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lleft[n][1]+x_lright[n-1][1])
						+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
						-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
						+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));	//Һ�崢��=����������+���洫��+������+��������-����
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
				/*�ҽ��������ƶ�һ�������壬n��n+1�ϲ�Ϊͬһ������*/
				else if (x_lleft[n][2]>=dx&&x_lleft[n][2]<2*dx)		
				{
					if (n==length)	//����������ҽ�������һ���������ƶ�����һ��������
					{
						/*�ҽ������������*/
                        T_lleft[n][2]=1/x_lleft[n][2]*(x_lleft[n][1]*T_lleft[n][1]+dt/(c_pl_left[n]*ro_l*A_l)
							*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]*(x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lleft[n][1]+x_lright[n-1][1])
							+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
							-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
							+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));	//ͬ�ϣ�ȫ��û�п���n+1�����壿����
						if (T_lleft[n][2]>=T_max)
						{
							T_lleft[n][2]=T_max;
						}
						if (T_lleft[n][2]<=T_min)
						{
							T_lleft[n][2]=T_min;
						}

						if (sort[1]!=4)	//��1�������岻Ϊ��-Һ����ֻ��Ϊ��
						{
							sort[1]=3;	//���Һ-��
                            x_lleft[1][2]=x_lleft[n][2]-dx;
                            x_lright[1][2]=0;
                            x_v[1][2]=dx-x_lleft[1][2]-x_lright[1][2];
                            v_lleft[1][2]=v_lleft[n][2];
                            v_lright[1][2]=0;
                            T_lleft[1][2]=T_lleft[n][2];
                            T_lright[1][2]=0;
						}
						else  //��1��������Ϊ��-Һ
						{
							sort[1]=8;	//���Һ-��-Һ
                            x_lleft[1][2]=x_lleft[n][2]-dx;
                            x_v[1][2]=dx-x_lleft[1][2]-x_lright[1][2];
                            v_lleft[1][2]=v_lleft[n][2];
                            T_lleft[1][2]=T_lleft[n][2];
						}
					}
					else  //�ҽ��治�����һ��������
					{
						/*�ҽ������������*/
                        T_lleft[n][2]=1/x_lleft[n][2]*(x_lleft[n][1]*T_lleft[n][1]+dt/(c_pl_left[n]*ro_l*A_l)
							*(v_lleft[n][1]*A_l*ro_l*c_pl_left[n]*(x_lright[n-1][1]*T_lleft[n][1]+x_lleft[n][1]*T_lright[n-1][1])/(x_lleft[n][1]+x_lright[n-1][1])
							+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]
							-2*A_l*(T_lleft[n][1]-T_lright[n-1][1])/(x_lleft[n][1]/lamt_left[n]+x_lright[n-1][1]/lamt_right[n-1])
							+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));	//ͬ���ϣ�ȫ��û�п���n+1�����壿����
						if (T_lleft[n][2]>=T_max)
						{
							T_lleft[n][2]=T_max;
						}
						if (T_lleft[n][2]<=T_min)
						{
							T_lleft[n][2]=T_min;
						}

						if (sort[n+1]!=4)	//������-Һ��ֻ������
						{
							sort[n+1]=3;	//���Һ-��
                            x_lleft[n+1][2]=x_lleft[n][2]-dx;
                            x_lright[n+1][2]=0;
                            x_v[n+1][2]=dx-x_lleft[n+1][2]-x_lright[n+1][2];
                            v_lleft[n+1][2]=v_lleft[n][2];
                            v_lright[n+1][2]=0;
                            T_lleft[n+1][2]=T_lleft[n][2];
                            T_lright[n+1][2]=0;
						}
						else  //��-Һ
						{
							sort[n+1]=8;	//���Һ-��-Һ
                            x_lleft[n+1][2]=x_lleft[n][2]-dx;
                            x_v[n+1][2]=dx-x_lleft[n+1][2]-x_lright[n+1][2];
                            v_lleft[n+1][2]=v_lleft[n][2];
                            T_lleft[n+1][2]=T_lleft[n][2];
						}
					}

					sort[n]=1;	//ԭ�ҽ����������Һ
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
				/*�ҽ��������ƶ�һ�������壬n��n-1�ϲ�Ϊͬһ������*/
				else		
				{
					x_lleft[n-1][2]=dx+x_lleft[n][2];
                    x_lright[n-1][2]=0;
                    x_v[n-1][2]=dx-x_lleft[n-1][2];
                    T_lright[n-1][2]=0;
                    v_lright[n-1][2]=0;
                    sort[n-1]=3;	//���Һ-��

					/*�ҽ������������*/
                    T_lleft[n-1][2]=1/(x_lleft[n-1][2]*c_pl_left[n-1])*(x_lleft[n][1]*T_lleft[n][1]*c_pl_left[n]+x_lright[n-1][1]*T_lright[n-1][1]*c_pl_right[n-1]
						+dt/(A_l*ro_l)*(v_lright[n-1][1]*A_l*ro_l*c_pl_right[n-1]*(x_lright[n-2][1]*T_lleft[n-1][1]+x_lleft[n-1][1]*T_lright[n-2][1])/(x_lleft[n-1][1]+x_lright[n-2][1])
						+h_lleft[n]*(T_w[n][1]-T_lleft[n][1])*pi*d*x_lleft[n][1]+h_lright[n-1]*(T_w[n-1][1]-T_lright[n-1][1])*pi*d*x_lright[n-1][1]
						-2*A_l*(T_lleft[n-1][1]-T_lright[n-2][1])/(x_lleft[n-1][1]/lamt_left[n-1]+x_lright[n-2][1]/lamt_right[n-2])
						+dl_right[n_liquid]*A_l*ro_l*c_pl_left[n]*T_lleft[n][1]/dt-dl_vl_right[n_liquid]*A_l*ro_l/dt*c_pl_left[n]*T_lleft[n][1]));	//������n-1��n������
					if (T_lleft[n-1][2]>=T_max)
					{
						T_lleft[n-1][2]=T_max;
					}
					if (T_lleft[n-1][2]<=T_min)
					{
						T_lleft[n-1][2]=T_min;
					}

					if (sort[n]==8)	//ԭ�ҽ��������ΪҺ-��-Һ
					{
						sort[n]=4;	//�����-Һ
                        T_lleft[n][2]=0;
                        v_lleft[n][2]=0;
                        x_lleft[n][2]=0;
                        x_v[n][2]=dx-x_lright[n][2];
					}
					else  //ԭ�ҽ��������ΪҺ-��
					{
						sort[n]=2;	//�����
                        T_lleft[n][2]=0;
                        T_lright[n][2]=0;
                        v_lleft[n][2]=0;
                        v_lright[n][2]=0;
                        x_lleft[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx;
					}
				}
				/*Һ��������浽�ҽ���������ϣ�������һ��Һ��*/
				n_liquid=n_liquid+1;
			}

			/*��������Χ���������ѭ��while (i<=nn)*/
			if (codenumber==1)
				{break;}
		}

		/******������ʼ��״̬�ķ������ۣ��ڶ����֣�Һ��Һ-����Һ-��-Һ����ʱ�������������Һ�����******/
		else if (sort[1]==1||sort[1]==3||sort[1]==8)
		{
			/*ȷ����Һ���������仯���������á�Һ���ܳ���*/

			/***********************�ȶ�n_liquid==1�����ۣ���1��Һ�����������ԭ��************************/
			n_liquid=1;
			jump=0;
			dm_left[n_liquid]=-Tr_v[n_liquid_total]/2*dt;	//���һ���������������һ��
			dl_left[n_liquid]=dm_left[n_liquid]/A_l/ro_l;	//��߽��ƶ�����
			dm_right[n_liquid]=-Tr_v[n_liquid]/2*dt;	//��1���������������һ��
			dl_right[n_liquid]=dm_right[n_liquid]/A_l/ro_l;	//�ұ߽��ƶ�����
        
			n=n_trans[2*n_liquid_total];	//��λ�����һ����Һ���棬����1��Һ���������
			G_force[n_liquid]=0;
			l[n_liquid]=0;
			C_ll=0;
			K=0;
        
			G_force[n_liquid]=G_force[n_liquid]+ro_l*g*c_gravity[n]*x_lright[n][1];	//��������������ۼ�
			l[n_liquid]=l[n_liquid]+x_lright[n][1];
			Re_lright[n]=d*ro_l*fabs(v_lright[n][1])/mu_right[n];
			if (Re_lright[n]==0)
				{C_l[n]=0;}
			else if (Re_lright[n]<=1180)
				{C_l[n]=16/Re_lright[n];}
			else 
				{C_l[n]=0.078*pow(Re_lright[n],-0.2);}
			if(Re_lright[n]>0&&(n%210>=101)&&(n%210<=210))	//��ͷ����������ʧ
				{K=1000/Re_lright[n]+0.1*(1+4/pow(2.3/25.4,0.3));
				C_l[n]=C_l[n]+K/4*2.3/110;}
			if (v_lright[n_trans[2*n_liquid-1]][1]>0)	//�ٶ���������������
				{C_l[n]=-C_l[n];}
			C_ll=C_ll+C_l[n]*x_lright[n][1];	//Һ��������ϵ��Ϊ������������ϵ���ļ�Ȩƽ��

			n=n+1;
			/*������浽����ԭ���ⲿ��Һ��*/
			while (n<=length)
			{
				G_force[n_liquid]=G_force[n_liquid]+ro_l*g*c_gravity[n]*x_lleft[n][1];	//��������
				l[n_liquid]=l[n_liquid]+x_lleft[n][1];	//Һ������
				Re_lleft[n]=d*ro_l*fabs(v_lleft[n][1])/mu_left[n];
				if (Re_lleft[n]==0)
					{C_l[n]=0;}
				else if (Re_lleft[n]<=1180)
					{C_l[n]=16/Re_lleft[n];}
				else 
					{C_l[n]=0.078*pow(Re_lleft[n],-0.2);}
				if(Re_lleft[n]>0&&(n%210>=101)&&(n%210<=210))	//��ͷ����������ʧ
					{K=1000/Re_lleft[n]+0.1*(1+4/pow(2.3/25.4,0.3));
					C_l[n]=C_l[n]+K/4*2.3/110;}
				if (v_lleft[n_trans[2*n_liquid-1]][1]>0)
					{C_l[n]=-C_l[n];}
				C_ll=C_ll+C_l[n]*x_lleft[n][1];
				n=n+1;
			}
			n=1;
			/*������ԭ�㵽�ҽ����ⲿ��Һ��*/
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
				if(Re_lleft[n]>0&&(n%210>=101)&&(n%210<=210))	//��ͷ����������ʧ
					{K=1000/Re_lleft[n]+0.1*(1+4/pow(2.3/25.4,0.3));
					C_l[n]=C_l[n]+K/4*2.3/110;}
				if (v_lleft[n_trans[2*n_liquid-1]][1]>0)
					{C_l[n]=-C_l[n];}
				C_ll=C_ll+C_l[n]*x_lleft[n][1];
				n=n+1;
			}

			C_ll=C_ll/l[n_liquid];	//֮ǰ��Ȩ�ۼӣ��˴�ƽ��

			/**�������̣����߶�������Һ�������A_l���õ��ҽ����������ٶ�**/
			v_lleft[n_trans[2*n_liquid-1]][2]=1/(l[n_liquid]+dl_left[n_liquid]+dl_right[n_liquid]-dl_vl_left[n_liquid]-dl_vl_right[n_liquid])
				*(v_lleft[n_trans[2*n_liquid-1]][1]*l[n_liquid]+dt/ro_l*(G_force[n_liquid]+P_v[n_trans[2*n_liquid_total]][1]-P_v[n_trans[2*n_liquid-1]][1]
				+2*C_ll*ro_l/d*l[n_liquid]*pow(v_lleft[n_trans[2*n_liquid-1]][1],2)));

			/*���������v_lֵ��������ͬһҺ���Ŀ�����*/
			n=n_trans[2*n_liquid_total];	//��λ�������
			while (n<=length)	//������浽����ԭ��
			{
				if (n!=n_trans[2*n_liquid_total])
				{
					v_lleft[n][2]=v_lleft[n_trans[2*n_liquid-1]][2]; 
					v_lright[n][2]=v_lleft[n_trans[2*n_liquid-1]][2];
				}
				else  //�����û����߲���Һ��
				{
					v_lright[n][2]=v_lleft[n_trans[2*n_liquid-1]][2];
				}
				n=n+1;
			}

			n=1;	//��λ����1��������
			while (n<=n_trans[2*n_liquid-1])	//������ԭ�㵽�ҽ���
			{
				if (n!=n_trans[2*n_liquid-1])
				{
					v_lleft[n][2]=v_lleft[n_trans[2*n_liquid-1]][2];
					v_lright[n][2]=v_lleft[n_trans[2*n_liquid-1]][2];
				}
				else  //�ҽ���û���ұ߲���Һ��
				{
					v_lleft[n][2]=v_lleft[n_trans[2*n_liquid-1]][2];
				}
				n=n+1;
			}

			/*�����ж��Ƿ���λ�õ�ƫ�ƣ�Һ���������յ������λ�øı䣩*/

			/*�����*/
			n=n_trans[2*n_liquid_total];
			x_lright[n][2]=x_lright[n][1]-v_lright[n][1]*dt+dl_left[n_liquid]-dl_vl_left[n_liquid];	//�����������µ�Һ�峤��
			if (sort[n]!=8)	//�����Һ-��-Һ������߲���Һ�峤�ȱ��ֲ���
				{x_lleft[n][2]=x_lleft[n][1];}
			/*��������ƶ�����һ�������壬����ѭ��while (i<=nn)*/
			if (x_lright[n][2]<=-dx||x_lright[n][2]>=2*dx)	
			{
				codenumber=2;
				break;
			}

			/*����������������ƫ��*/
			if (x_lright[n][2]>0&&x_lright[n][2]<dx)			
			{
				/*�����Һ�����������*/
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
				else  //�������������һ�������壬��n+1=1
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
			/*����������ƶ�һ�������壬��n��n+1�ϲ�Ϊһ��������*/
			else if (x_lright[n][2]<=0)							
			{
				/*�����ԭ�������һ�������壬�ƶ�����1��������*/
				if (n==length)
				{
					x_lright[1][2]=dx+x_lright[n][2];
                    x_lleft[1][2]=0;
                    x_v[1][2]=dx-x_lright[1][2]-x_lleft[1][2];
                    T_lleft[1][2]=0;
                    v_lleft[1][2]=0;
                    sort[1]=4;	//��1������������-Һ
                    jump=1;
                    cross_positive=1;	//Һ�������������ԭ�㣬������ż�1

					/*�����Һ�����������*/
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

					if (sort[n]!=8)	//����������ԭ��Ϊ��-Һ
					{
						T_lleft[n][2]=0;
                        T_lright[n][2]=0;
                        v_lleft[n][2]=0;
                        v_lright[n][2]=0;
                        x_lleft[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx;
                        sort[n]=2;	//�����ƶ�������
					}
					else  //����������ԭ��ΪҺ-��-Һ
					{
						sort[n]=3;	//�����ƶ�����Һ-��
                        T_lright[n][2]=0;
                        v_lright[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx-x_lleft[n][2];
					}
				}
				/*�����ԭ���ڵ�����2�������壬�ƶ������1��������*/
				else if (n==(length-1))
				{
					x_lright[n+1][2]=dx+x_lright[n][2];
                    x_lleft[n+1][2]=0;
                    x_v[n+1][2]=dx-x_lright[n+1][2]-x_lleft[n+1][2];
                    T_lleft[n+1][2]=0;
                    v_lleft[n+1][2]=0;
                    sort[n+1]=4;	//���һ������������-Һ
                    jump=1;

					/*�����Һ�����������*/
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

					if (sort[n]!=8)	//����������ԭ��Ϊ��-Һ
					{
						T_lleft[n][2]=0;
                        T_lright[n][2]=0;
                        v_lleft[n][2]=0;
                        v_lright[n][2]=0;
                        x_lleft[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx;
                        sort[n]=2;	//���ڱ����
					}
					else//����������ԭ��ΪҺ-��-Һ
					{
						sort[n]=3;	//���ڱ��Һ-��
                        T_lright[n][2]=0;
                        v_lright[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx-x_lleft[n][2];
					}
				}
				/*�����ԭ���ڵ�����3����������߸���*/
				else
				{
					x_lright[n+1][2]=dx+x_lright[n][2];
                    x_lleft[n+1][2]=0;
                    x_v[n+1][2]=dx-x_lright[n+1][2]-x_lleft[n+1][2];
                    T_lleft[n+1][2]=0;
                    v_lleft[n+1][2]=0;
                    sort[n+1]=4;	//n+1������ڵ����������壬��-Һ
                    jump=1;

					/*�����Һ�����������*/
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

					if (sort[n]!=8)	//ԭ��Ϊ��-Һ
					{
						T_lleft[n][2]=0;
                        T_lright[n][2]=0;
                        v_lleft[n][2]=0;
                        v_lright[n][2]=0;
                        x_lleft[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx;
                        sort[n]=2;	//�ֱ����
					}
					else  //ԭ��ΪҺ-��-Һ
					{
						sort[n]=3;	//����ΪҺ-��
                        T_lright[n][2]=0;
                        v_lright[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx-x_lleft[n][2];
					}
				}
			}
			/*����������ƶ�һ�������壬��n��n-1�ϲ�Ϊһ��������*/
			else
			{
				/*�����Һ�����������*/
				/*�����ԭ�������һ�������壬�ƶ���������2��������*/
				if (n==length)
				{
					T_lright[n][2]=1/x_lright[n][2]*(x_lright[n][1]*T_lright[n][1]+dt/(c_pl_right[n]*ro_l*A_l)
						*(-v_lright[n][1]*A_l*ro_l*c_pl_right[n]*(x_lleft[1][1]*T_lright[n][1]+x_lright[n][1]*T_lleft[1][1])/(x_lright[n][1]+x_lleft[1][1])
						+h_lright[n]*(T_w[n][1]-T_lright[n][1])*pi*d*x_lright[n][1]
						-2*A_l*(T_lright[n][1]-T_lleft[1][1])/(x_lright[n][1]/lamt_right[n]+x_lleft[1][1]/lamt_left[1])
						+dl_left[n_liquid]*A_l*ro_l*c_pl_right[n]*T_lright[n][1]/dt-dl_vl_left[n_liquid]*A_l*ro_l/dt*c_pl_right[n]*T_lright[n][1]));
				}
				/*�����ԭ���ڵ�����2����������߸���*/
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

				if (sort[n-1]!=3)	//ԭ�������߿�����Ϊ��
				{
					x_lright[n-1][2]=x_lright[n][2]-dx; 
                    x_lleft[n-1][2]=0;
                    x_v[n-1][2]=dx-x_lright[n-1][2]-x_lleft[n-1][2];
                    v_lright[n-1][2]=v_lright[n][2];
                    v_lleft[n-1][2]=0;
                    T_lright[n-1][2]=T_lright[n][2];
                    T_lleft[n-1][2]=0;
                    sort[n-1]=4;	//�ֱ�Ϊ��-Һ
                          
                    x_lright[n][2]=dx;
                    x_lleft[n][2]=dx;
                    x_v[n][2]=0;
                    v_lleft[n][2]=v_lright[n][2];
                    T_lleft[n][2]=T_lright[n][2];
                    ro_v[n][2]=0;
                    T_v[n][2]=0;
					T_sat[n][2] = 0;
                    P_v[n][2]=0;
                    sort[n]=1;	//ԭ������ΪҺ
				}
				else //ԭ�������߿�����ΪҺ-��
				{
					x_lright[n-1][2]=x_lright[n][2]-dx; 
                    x_v[n-1][2]=dx-x_lright[n-1][2]-x_lleft[n-1][2];
                    v_lright[n-1][2]=v_lright[n][2];
                    T_lright[n-1][2]=T_lright[n][2];
                    sort[n-1]=8;	//�ֱ�ΪҺ-��-Һ

					x_lright[n][2]=dx;
                    x_lleft[n][2]=dx;
                    x_v[n][2]=0;
                    v_lleft[n][2]=v_lright[n][2];
                    T_lleft[n][2]=T_lright[n][2];
                    ro_v[n][2]=0;
                    T_v[n][2]=0;
					T_sat[n][2] = 0;
                    P_v[n][2]=0;
                    sort[n]=1;	//ԭ������ΪҺ
				}

			}

			/*Һ�����м䲿��*/
			if (jump==1)	//����������Ƴ�������
			{
				if (n_trans[2*n_liquid_total]==length)	//�����ԭ�������һ�������壬��������1��������
				{
					n=2;
					while (n<=n_trans[2*n_liquid-1]-1)	//�ӵ�2�������嵽�ҽ�����߿�����
					{
						/*Һ��������������̣�����ͬ����dx*/
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
				else if (n_trans[2*n_liquid_total]==(length-1))	//�����ԭ���ڵ�����2�������壬���ƶ������һ��������
				{
					n=1;
					while (n<=n_trans[2*n_liquid-1]-1)	//�ӵ�1�������嵽�ҽ�����߿�����
					{
						if (n==1)	//��1�������嵥������
						{
							/*Һ���������������*/
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
						else  //�ӵ�2�������嵽�ҽ�����߿�����
						{
							/*Һ���������������*/
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
				else  //�����ԭ���ڵ�����3����������߸��󣬽����Һ����ԭ��Ϊ���Ϊ����
				{
					n=n_trans[2*n_liquid_total]+2;	//�µ�������ұ߿����壬�����һ��������
					while (n<=length)
					{
						if (n==length)	//���һ�������嵥������
						{
							/*Һ������������*/
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
						else  //�µ�������ұ߿����壬��������2��������
						{
							/*Һ������������*/
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
					n=1;	//�ӵ�1�������嵽�ҽ���
					if (n_trans[1]>=2)	//�ҽ����ڵ�2����������߸��ң������1������������ҽ��棬�ں������
					{
						while (n<=n_trans[1]-1)	//�ӵ�1�������嵽�ҽ�����߿�����
						{
							if (n==1)	//��1�������嵥������
							{
								/*Һ������������*/
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
							else //�ӵ�2�������嵽�ҽ�����߿�����
							{
								/*Һ������������*/
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
			else  //�������ԭ��������������Ƴ�������
			{
				if (n_trans[2*n_liquid_total]==length)	//����������һ��������
				{
					n=1;
					while (n<=n_trans[2*n_liquid-1]-1)	//�ӵ�1�������嵽�ҽ�����߿�����
					{
						if (n==1)	//��1�������嵥������
						{
							/*Һ������������*/
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
						else  //�ӵ�2�������嵽�ҽ�����߿�����
						{
							/*Һ������������*/
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
				else  //����治�����һ�������壬��Һ����ԭ��Ϊ��ֶ�����
				{
					n=n_trans[2*n_liquid_total]+1;
					while (n<=length)	//��������ұ߿����嵽���һ��������
					{
						if (n==length)	//���һ�������嵥������
						{
							/*Һ������������*/
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
						else  //��������ұ߿����嵽������2��������
						{
							/*Һ������������*/
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
					n=1;	//�ӵ�1�������嵽�ҽ�����߿�����
					if (n_trans[1]>=2)	//�ұ߽��ڵ�2����������߸��ұߣ������1������������ұ߽磬�ں������
					{
						while (n<=n_trans[1]-1)
						{
							if (n==1)	//��1�������嵥������
							{
								/*Һ������������*/
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
							else  //�ӵ�2�������嵽�ҽ�����߿�����
							{
								/*Һ������������*/
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

			jump=0;	//���ý����ƶ����

			/*�ҽ���*/
			n=n_trans[2*n_liquid-1];
			x_lleft[n][2]=x_lleft[n][1]+v_lleft[n][1]*dt+dl_right[n_liquid]-dl_vl_right[n_liquid];	//�ҽ���������µ�Һ�峤��
			if (sort[n]!=8)	//�ҽ��治��Һ-��-Һ������Һ-��
				{x_lright[n][2]=x_lright[n][1];}
			/*���ҽ����ƶ�����һ�������壬����ѭ��while (i<=nn)*/
			if (x_lleft[n][2]<=-dx||x_lleft[n][2]>=2*dx)
			{
				codenumber=2;
				break;
			}
			/*�ҽ���������ƫ��*/
			if (x_lleft[n][2]>0&&x_lleft[n][2]<dx)
			{
				if (n==1)//�ҽ����ڵ�1��������
				{
					/*�ҽ���Һ������������*/
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
				else//�ҽ����ڵ�2������������
				{
					/*�ҽ���Һ������������*/
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
			/*�ҽ��������ƶ�һ�������壬n��n+1�ϲ�����*/
			else if (x_lleft[n][2]>=dx&&x_lleft[n][2]<2*dx)		
			{
				/*�ҽ���Һ������������*/
				if (n==1)//�ҽ���ԭ���ڵ�1��������
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
				else//�ҽ���ԭ���ڵ�2������������
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

				if (sort[n+1]!=4)//ԭ��Ϊ��
				{
					sort[n+1]=3;	//���ڱ��Һ-��
                    x_lleft[n+1][2]=x_lleft[n][2]-dx;
                    x_lright[n+1][2]=0;
                    x_v[n+1][2]=dx-x_lleft[n+1][2]-x_lright[n+1][2];
                    v_lleft[n+1][2]=v_lleft[n][2];
                    v_lright[n+1][2]=0;
                    T_lleft[n+1][2]=T_lleft[n][2];
                    T_lright[n+1][2]=0;
				}
				else//ԭ��Ϊ��-Һ
				{
					sort[n+1]=8;	//���ڱ��Һ-��-Һ
					x_lleft[n+1][2]=x_lleft[n][2]-dx;
					x_v[n+1][2]=dx-x_lleft[n+1][2]-x_lright[n+1][2];
					v_lleft[n+1][2]=v_lleft[n][2];
					T_lleft[n+1][2]=T_lleft[n][2];
				}
				sort[n]=1;	//ԭ�ҽ�����Һ
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
			/*�ҽ�������ƫ��һ�������壬n��n-1�ϲ�����*/
			else												
			{
				if (n==1)	//ԭ�ҽ����ڵ�1�������壬���ƶ������һ��������
				{
					x_lleft[length][2]=dx+x_lleft[n][2];
                    x_lright[length][2]=0;
                    x_v[length][2]=dx-x_lleft[length][2]-x_lright[length][2];
                    v_lright[length][2]=0;
                    T_lright[length][2]=0;
                    sort[length]=3;	//Һ-��
				}
				else//ԭ�ҽ����ڵ�2����������߸���
				{
					x_lleft[n-1][2]=dx+x_lleft[n][2];
                    x_lright[n-1][2]=0;
                    x_v[n-1][2]=dx-x_lleft[n-1][2]-x_lright[n-1][2];
                    v_lright[n-1][2]=0;
                    T_lright[n-1][2]=0;
                    sort[n-1]=3;	//Һ-��
				}
				/*�ҽ���Һ������������*/
				if (n==1)	//ԭ�ҽ����ڵ�1�������壬���ƶ������һ��������
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
				else if (n==2)	//ԭ�ҽ����ڵ�2�������壬���ƶ�����1��������
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
				else//ԭ�ҽ����ڵ�3����������߸��ң����ƶ�����2������������
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
				if (sort[n]==8)	//ԭ�ҽ���ΪҺ-��-Һ
				{
					sort[n]=4;	//�ֱ�Ϊ��-Һ
                    x_lleft[n][2]=0;
                    x_v[n][2]=dx-x_lleft[n][2]-x_lright[n][2];
                    T_lleft[n][2]=0;
                    v_lleft[n][2]=0;
				}
				else//ԭ�ҽ���ΪҺ-��
				{
					sort[n]=2;	//�ֱ�Ϊ��
                    T_lleft[n][2]=0;
                    T_lright[n][2]=0;
                    v_lleft[n][2]=0;
                    v_lright[n][2]=0;
                    x_lleft[n][2]=0;
                    x_lright[n][2]=0;
                    x_v[n][2]=dx;
				}
			}

			/*%%%%%%%%%%%%%%%%%%%%%��n_liquid!=1ʱ������%%%%%%%%%%%%%%%%%%%%%%%%%*/
			n_liquid=2;
			while (n_liquid<=n_liquid_total)
			{
				dm_left[n_liquid]=-Tr_v[n_liquid-1]/2*dt;	//������ݷ�һ������
				dl_left[n_liquid]=dm_left[n_liquid]/A_l/ro_l;	//������ƶ�����
				dm_right[n_liquid]=-Tr_v[n_liquid]/2*dt;	//�ұ����ݷ�һ������
				dl_right[n_liquid]=dm_right[n_liquid]/A_l/ro_l;    //�ҽ����ƶ�����
				n=n_trans[2*n_liquid-2];	//�����
				G_force[n_liquid]=0;
				l[n_liquid]=0;
				C_ll=0;
				K=0;
             
				G_force[n_liquid]=G_force[n_liquid]+ro_l*g*c_gravity[n]*x_lright[n][1];	//�������������ۼ�
				l[n_liquid]=l[n_liquid]+x_lright[n][1];	//��������Һ�峤����x_lright[n][1]��ʾ
				Re_lright[n]=d*ro_l*fabs(v_lright[n][1])/mu_right[n];
				if (Re_lright[n]==0)
					{C_l[n]=0;}
				else if (Re_lright[n]<=1180)
					{C_l[n]=16/Re_lright[n];}
				else
					{C_l[n]=0.078*pow(Re_lright[n],-0.2);}
				if(Re_lright[n]>0&&(n%210>=101)&&(n%210<=210))	//��ͷ����������ʧ
					{K=1000/Re_lright[n]+0.1*(1+4/pow(2.3/25.4,0.3));
					C_l[n]=C_l[n]+K/4*2.3/110;}
				if (v_lleft[n][1]>0)
					{C_l[n]=-C_l[n];}
				C_ll=C_ll+C_l[n]*x_lright[n][1];	//����ϵ���ȼ�Ȩ�ۼ�

				n=n+1;
				while (n<=n_trans[2*n_liquid-1])	//��������ұ߿����嵽�ҽ�����߿�����
				{
					G_force[n_liquid]=G_force[n_liquid]+ro_l*g*c_gravity[n]*x_lleft[n][1];	//�����ۼ�
					l[n_liquid]=l[n_liquid]+x_lleft[n][1];	//Һ�峤���ۼ�
					Re_lleft[n]=d*ro_l*fabs(v_lleft[n][1])/mu_left[n];
					if (Re_lleft[n]==0)
						{C_l[n]=0;}
					else if (Re_lleft[n]<=1180)
						{C_l[n]=16/Re_lleft[n];}
					else
						{C_l[n]=0.078*pow(Re_lleft[n],-0.2);}
					if(Re_lright[n]>0&&(n%210>=101)&&(n%210<=210))	//��ͷ����������ʧ
					{K=1000/Re_lleft[n]+0.1*(1+4/pow(2.3/25.4,0.3));
					C_l[n]=C_l[n]+K/4*2.3/110;}
					if (v_lleft[n][1]>0)
						{C_l[n]=-C_l[n];}
					C_ll=C_ll+C_l[n]*x_lleft[n][1];	//����ϵ����Ȩ�ۼ�
					n=n+1;
				}
				C_ll=C_ll/l[n_liquid];	//ǰ���Ȩ�ۼӣ��˴�ƽ��

				/*Һ���Ķ������̣��õ��ҽ����ٶ�*/
				v_lleft[n_trans[2*n_liquid-1]][2]=1/(l[n_liquid]+dl_left[n_liquid]+dl_right[n_liquid]-dl_vl_left[n_liquid]-dl_vl_right[n_liquid])
					*(v_lleft[n_trans[2*n_liquid-1]][1]*l[n_liquid]+dt/ro_l*(G_force[n_liquid]+P_v[n_trans[2*n_liquid-2]][1]-P_v[n_trans[2*n_liquid-1]][1]
					+2*C_ll*ro_l/d*l[n_liquid]*pow(v_lleft[n_trans[2*n_liquid-1]][1],2)));

				/*���������v_lֵ��������ͬһҺ���Ŀ�����*/
				n=n_trans[2*n_liquid-2];	//�����
				while (n<=n_trans[2*n_liquid-1]-1)	//������浽�ҽ������
				{
					if (n==n_trans[2*n_liquid-2])	//�����
					{
						v_lright[n][2]=v_lleft[n_trans[2*n_liquid-1]][2];
					}
					else//���������
					{
						v_lleft[n][2]=v_lleft[n_trans[2*n_liquid-1]][2]; 
						v_lright[n][2]=v_lleft[n_trans[2*n_liquid-1]][2];
					}
					n=n+1;
				}

				/**�����ж��Ƿ���λ�õ�ƫ�ƣ�Һ���������յ������λ�øı䣩**/
				/*�����*/
				n=n_trans[2*n_liquid-2];
				x_lright[n][2]=x_lright[n][1]-v_lright[n][1]*dt+dl_left[n_liquid]-dl_vl_left[n_liquid];	//������µ�Һ�峤��
				if (sort[n]!=8)	//�����Ϊ��-Һ
					{x_lleft[n][2]=x_lleft[n][1];}
				/*��������ƶ�����һ�������壬����ѭ��while (n_liquid<=n_liquid_total)��֮��Ҫ�������ѭ��*/
				if (x_lright[n][2]<=-dx||x_lright[n][2]>=2*dx)
				{
					codenumber=2;
					break;
				}
				/*�������ƫ�ƣ�����ԭ������*/
				if (x_lright[n][2]>0&&x_lright[n][2]<dx)
				{
					/*%%Һ�����������%%%*/
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
				/*����������ƶ�һ�������壬��n��n+1�ϲ�����*/
				else if (x_lright[n][2]<=0)								
				{
					x_lright[n+1][2]=dx+x_lright[n][2];
                    x_lleft[n+1][2]=0;
                    x_v[n+1][2]=dx-x_lright[n+1][2]-x_lleft[n+1][2];
                    T_lleft[n+1][2]=0;
                    v_lleft[n+1][2]=0;
                    sort[n+1]=4;	//��-Һ
                    jump=1;	//��������Ʊ��

					/*�����Һ�����������*/
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
					if (sort[n]!=8)	//ԭ�����Ϊ��-Һ
					{
						T_lleft[n][2]=0;
                        T_lright[n][2]=0;
                        v_lleft[n][2]=0;
                        v_lright[n][2]=0;
                        x_lleft[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx;
                        sort[n]=2;	//�ֱ����
					}
					else//ԭ�����ΪҺ-��-Һ
					{
						sort[n]=3;	//�ֱ�ΪҺ-��
                        T_lright[n][2]=0;
                        v_lright[n][2]=0;
                        x_lright[n][2]=0;
                        x_v[n][2]=dx-x_lleft[n][2];
					}
				}
				/*����������ƶ�һ�������壬��n��n-1�ϲ�����*/
				else															
				{
					/*�����Һ�����������*/
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

					if (sort[n-1]!=3)	//ԭ��������Ϊ��
					{
						x_lright[n-1][2]=x_lright[n][2]-dx; 
                        x_lleft[n-1][2]=0;
                        x_v[n-1][2]=dx-x_lright[n-1][2]-x_lleft[n-1][2];
                        v_lright[n-1][2]=v_lright[n][2];
                        v_lleft[n-1][2]=0;
                        T_lright[n-1][2]=T_lright[n][2];
                        T_lleft[n-1][2]=0;
                        sort[n-1]=4;	//�ֱ�Ϊ��-Һ
                          
                        x_lright[n][2]=dx;
                        x_lleft[n][2]=dx;
                        x_v[n][2]=0;
                        v_lleft[n][2]=v_lright[n][2];
                        T_lleft[n][2]=T_lright[n][2];
                        ro_v[n][2]=0;
                        T_v[n][2]=0;
						T_sat[n][2] = 0;
                        P_v[n][2]=0;
                        sort[n]=1;	//ԭ������ΪҺ
					}
					else//ԭ�����ΪҺ-��
					{
						x_lright[n-1][2]=x_lright[n][2]-dx; 
                        x_v[n-1][2]=dx-x_lright[n-1][2]-x_lleft[n-1][2];
                        v_lright[n-1][2]=v_lright[n][2];
                        T_lright[n-1][2]=T_lright[n][2];
                        sort[n-1]=8;	//�ֱ�ΪҺ-��-Һ
                          
                        x_lright[n][2]=dx;
                        x_lleft[n][2]=dx;
                        x_v[n][2]=0;
                        v_lleft[n][2]=v_lright[n][2];
                        T_lleft[n][2]=T_lright[n][2];
                        ro_v[n][2]=0;
                        T_v[n][2]=0;
						T_sat[n][2] = 0;
						P_v[n][2]=0;
						sort[n]=1;	//ԭ������ΪҺ
					}
				}

				/*Һ�����м䲿��*/
				if (jump==1)//���������
				{
					n=n_trans[2*n_liquid-2]+2;	//��λ����������ұ߿�����
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
				else//��������ƻ�����ԭ������
				{
					n=n_trans[2*n_liquid-2]+1;	//��λ��ԭ������ұ߿�����
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
				jump=0;	//����������ƶ����

				/*�ҽ���*/
				n=n_trans[2*n_liquid-1];
				x_lleft[n][2]=x_lleft[n][1]+v_lleft[n][1]*dt+dl_right[n_liquid]-dl_vl_right[n_liquid];	//�ҽ���������µ�Һ�峤��
				if (sort[n]!=8)	//�ҽ���ΪҺ-��
					{x_lright[n][2]=x_lright[n][1];}
				/*���ҽ����ƶ�����һ�������壬����ѭ��while (n_liquid<=n_liquid_total)��֮��Ҫ�������ѭ��*/
				if (x_lleft[n][2]<=-dx||x_lleft[n][2]>=2*dx)
				{
					codenumber=2;
					break;
				}
				/*�ҽ���û���Ƴ�ԭ������*/
				if (x_lleft[n][2]>0&&x_lleft[n][2]<dx)
				{
					/*�ҽ���Һ������������*/
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
				/*�ҽ��������ƶ�һ�������壬n��n+1�ϲ�����*/
				else if (x_lleft[n][2]>=dx&&x_lleft[n][2]<2*dx)
				{
					/*�ҽ���Һ������������*/
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

					if (sort[n+1]!=4)	//ԭ�ҽ����ұ߿�����Ϊ��
					{
						sort[n+1]=3;	//�ֱ��Һ-��
                        x_lleft[n+1][2]=x_lleft[n][2]-dx;
                        x_lright[n+1][2]=0;
                        x_v[n+1][2]=dx-x_lleft[n+1][2]-x_lright[n+1][2];
                        v_lleft[n+1][2]=v_lleft[n][2];
                        v_lright[n+1][2]=0;
                        T_lleft[n+1][2]=T_lleft[n][2];
                        T_lright[n+1][2]=0;
					}
					else//ԭ�ҽ����ұ߿�����Ϊ��-Һ
					{
						sort[n+1]=8;	//�ֱ��Һ-��-Һ
                        x_lleft[n+1][2]=x_lleft[n][2]-dx;
                        x_v[n+1][2]=dx-x_lleft[n+1][2]-x_lright[n+1][2];
                        v_lleft[n+1][2]=v_lleft[n][2];
                        T_lleft[n+1][2]=T_lleft[n][2];
					}
					sort[n]=1;	//ԭ�ҽ�����Һ
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
				/*�ҽ��������ƶ�һ�������壬n��n-1�ϲ�����*/
				else													
				{
					x_lleft[n-1][2]=dx+x_lleft[n][2];
                    x_lright[n-1][2]=0;
                    x_v[n-1][2]=dx-x_lleft[n-1][2]-x_lright[n-1][2];
                    v_lright[n-1][2]=0;
                    T_lright[n-1][2]=0;
                    sort[n-1]=3;	//���ҽ���Һ-��
					
					/*�ҽ���Һ�����������*/
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

					if (sort[n]==8)	//ԭ�ҽ���ΪҺ-��-Һ
					{
						sort[n]=4;	//�ֱ�Ϊ��-Һ
                        x_lleft[n][2]=0;
                        x_v[n][2]=dx-x_lleft[n][2]-x_lright[n][2];
                        T_lleft[n][2]=0;
                        v_lleft[n][2]=0;
					}
					else//ԭ�ҽ���ΪҺ-��
					{
						sort[n]=2;	//�ֱ�Ϊ��
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

			/*�����ҽ��泬����Χ���������ѭ��*/
			if (codenumber==2)
				{break;}
		}

		/******������ʼ��״̬�ķ������ۣ��������֣������쳣��������-Һ-��֮��*********/
		else
		{
			fid0=fopen("state.txt","a+");	//����a+ģʽ��һ������state.txt���ļ�������state.txt���ļ�ָ���fid0,�쳣�򷵻�NULL
			fprintf(fid0,"Strange,sort[1] impossible, n=%d, i=%d",n,i);	//����Strange,sort[1] impossible���Լ�n��iд����fid0ָ�����ļ�
			fclose(fid0);	//�ر���fid0ָ�����ļ�,���ز��������0��EOF
		}

		/*%%%%%%%%%%%%%%%%%%%������ͨ��x=0���µ�������ŵı仯����Һ����Ų���%%%%%%%%%%%%%%%%%%%%%%%%*/
		double P_vb_tem[100][3];	//������ű仯����������ֵ������
		double T_vb_tem[100][3];	//������ű仯����������ֵ������
		double T_sat_vb_tem[100][3];	//������ű仯����������ֵ������
		double ro_vb_tem[100][3];	//������ű仯����������ֵ������
		double V_bubble_tem[100][3];	//������ű仯����������ֵ������
		double m_bubble_tem[100][3];	//������ű仯����������ֵ������
		double ro_vb[100][3];	//���������ܶ�

		if (cross_negative==1)	//������ż�1
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
			P_vb_tem[n_bubble_total][1]=P_vb[1][1];	//��1�����ݣ�������һ��
			T_vb_tem[n_bubble_total][1]=T_vb[1][1];
			T_sat_vb_tem[n_bubble_total][1]=T_sat_vb[1][1];
			ro_vb_tem[n_bubble_total][1]=ro_vb[1][1];
			V_bubble_tem[n_bubble_total][1]=V_bubble[1][1];
			m_bubble_tem[n_bubble_total][1]=m_bubble[1][1];

			n_bubble=1;
			while (n_bubble<=n_bubble_total)	//��ű仯���״̬����
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

		if (cross_positive==1)	//������ż�1
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
			P_vb_tem[1][1]=P_vb[n_bubble_total][1];	//���һ�����ݱ�ɵ�1������
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


		/*%%%%�Գ���sort(n)==9(�������ҽ�������ͬһ�����壺ͨ������1113111��1114111)%%%*/

		int n_sort9;	//sort[n]==9�Ŀ��������
		int n_bubble_tem;
		n=1;
		if ((sort[1]==3||sort[1]==4)&&sort[length]*sort[2]==1)	//�������ҽ��涼�ڵ�1��������
			{sort[1]=9;}
		n=length;
		if ((sort[length]==3||sort[length]==4)&&sort[length-1]*sort[1]==1)	//�������ҽ��涼�����һ��������
			{sort[length]=9;}
		n=2;
		while (n<=length-1)
		{
			if ((sort[n]==3||sort[n]==4)&&sort[n-1]*sort[n+1]==1)
				{sort[n]=9;}
			n=n+1;
		}
        
		/*��������������sort[n]==9�ĵ����������*/
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
			fid0=fopen("state.txt","a+");	//����a+ģʽ��һ������state.txt���ļ�������state.txt���ļ�ָ���fid0,�쳣�򷵻�NULL
			fprintf(fid0,"Special, Codenumber=%d, n=%d, i=%d, n_sort9=%d",codenumber,n,i,n_sort9);	//������д����fid0ָ�����ļ�����ʱcodenumber==������
			fclose(fid0);	//�ر���fid0ָ�����ļ�,���ز��������0��EOF
		}
		codenumber=0;

		/*����sort[n]==9�����������ʽ�仯*/
		if (n_sort9!=0)
		{
			if (sort[1]==1)	//��1��������ΪҺ
			{
				n_bubble=1;	//��ʼΪ1����sort[n]==9һ������
				n=2;
				transfer_point_number=0;
				while (n<=length)
				{
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=transfer_point_number+1;
						if ((transfer_point_number%2)==0)	//����2������0
							{n_bubble=n_bubble+1;}
					}
					if (sort[n]==8)
					{
						transfer_point_number=transfer_point_number+2;
						n_bubble=n_bubble+1;
					}
					if (sort[n]==9)	//��sort[n]=9�����sort[n]=1��������ż�С1
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
			if (sort[1]==2)	//��1��������Ϊ��
			{
				n_bubble=1;	//��ǰ��sort[n]==9һ����������
				n=2;
				transfer_point_number=0;
				while (n<=length)
				{
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=transfer_point_number+1;
						if ((transfer_point_number%2)!=0)	//��������2
							{n_bubble=n_bubble+1;}
					}
					if (sort[n]==8)
					{
						transfer_point_number=transfer_point_number+2;
						n_bubble=n_bubble+1;
					}
					if (sort[n]==9)	//�����sort[n]==1��������ż�С1
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
						while (n_bubble_tem<=n_bubble_total-1)	//sort[n]==9֮���������ż�С1
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
			if (sort[1]==3)	//��1��������ΪҺ-��
			{
				n_bubble=1;
				n=2;
				transfer_point_number=1;
				while (n<=length)
				{
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=transfer_point_number+1;
						if ((transfer_point_number%2)==0)	//����2
							{n_bubble=n_bubble+1;}
					}
					if (sort[n]==8)
					{
						transfer_point_number=transfer_point_number+2;
						n_bubble=n_bubble+1;
					}
					if (sort[n]==9)//�����sort[n]==1��������ż�С1
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
			if (sort[1]==4)	//��1��������Ϊ��-Һ
			{
				n_bubble=2;
				n=2;
				transfer_point_number=1;
				while (n<=length)
				{
					if (sort[n]==3||sort[n]==4)
					{
						transfer_point_number=transfer_point_number+1;
						if ((transfer_point_number%2)!=0)	//��������2����Ϊ����
							{n_bubble=n_bubble+1;}
					}
					if (sort[n]==8)
					{
						transfer_point_number=transfer_point_number+2;
						n_bubble=n_bubble+1;
					}
					if (sort[n]==9)	//�����sort[n]==1��������ż�С1
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
			if (sort[1]==9)	//��1�����������һ��������
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
				sort[n]=1;	//�����sort[n]==1��������ż�С1

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
				transfer_point_number=2;	//����
				n_bubble=2;	//���࣬Ϊ�β���1������
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
					if (sort[n]==9)	//�����sort[n]==1��������ż�С1
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

		/*%%%%%%Һ���ϲ�λ�õ�ȷ�ϣ��·ֲ�������ֵ�����¶���%%%%%%%*/
		if (sort[1]==1)	//��1��������ΪҺ
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
				if (sort[n]==8)	//Һ���ڴ˿��ܺϲ�
				{
					x_v[n][2]=dx-x_lleft[n][2]-x_lright[n][2];
					if (x_v[n][2]<l_disappear)	//Һ���ϲ�
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
						sort[n]=1;	//�ϲ�����1���˺��������ż�С1
                   
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
					else  //Һ�����ϲ�
					{
						transfer_point_number=transfer_point_number+2;
						if ((transfer_point_number%2)==0)
							{n_bubble=n_bubble+1;}
					}
				}
				n=n+1;
			}
		}
		if (sort[1]==2)	//��1��������Ϊ��
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
						sort[n]=1;	//Һ���ϲ����˺�������ż�С1
                   
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
					else  //Һ�����ϲ�
					{
						transfer_point_number=transfer_point_number+2;
						if ((transfer_point_number%2)!=0)
							{n_bubble=n_bubble+1;}
					}
				}
				n=n+1;
			}
		}
		if (sort[1]==3)	//��1��������ΪҺ-��
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
					if (x_v[n][2]<l_disappear)	//Һ���ϲ�
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
					else  //Һ�����ϲ�
					{
						transfer_point_number=transfer_point_number+2;
						if ((transfer_point_number%2)==0)
							{n_bubble=n_bubble+1;}
					}
				}
				n=n+1;
			}

		}
		if (sort[1]==4)	//��1��������Ϊ��-Һ
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
					if (x_v[n][2]<l_disappear)	//Һ���ϲ�
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
					else  //Һ�����ϲ�
					{
						transfer_point_number=transfer_point_number+2;
						if ((transfer_point_number%2)!=0)
							{n_bubble=n_bubble+1;}
					}
				}
				n=n+1;
			}

		}
		if (sort[1]==8)	//��1��������ΪҺ-��-Һ
		{
			n_bubble=1;
			x_v[1][2]=dx-x_lleft[1][2]-x_lright[1][2];
			if (x_v[1][2]<l_disappear)	//Һ���ϲ�
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
			else  //Һ�����ϲ�
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
					if (x_v[n][2]<l_disappear)	//Һ���ϲ�
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
					else  //Һ�����ϲ�
					{
						transfer_point_number=transfer_point_number+2;
						if ((transfer_point_number%2)==0)
							{n_bubble=n_bubble+1;}
					}
				}
				n=n+1;
			}
		}
		/*����ȷ����Һ����λ�ú�����*/
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
				n_trans[transfer_point_number]=n;                  //ʹ��n_trans[]����¼ÿ���������λ��
			}
			if (sort[n]==8)
			{
				transfer_point_number=transfer_point_number+1;
				n_trans[transfer_point_number]=n;                  
				transfer_point_number=transfer_point_number+1;
				n_trans[transfer_point_number]=n;                  //n_trans{]����¼ÿ���������λ��
			}
			n=n+1;
		}

		n_trans_total=transfer_point_number;	//��¼��Һ��������
		n_bubble_total=n_trans_total/2;	//���ݵ�����
		n_liquid_total=n_trans_total/2;  //Һ��������

		/*%%%%%%%%����Һ���ϲ�֮���Һ���ٶ�%%%%%%%%*/
		double momentum;	//��λ������m2/s
		double l_liquid;	//Һ�����ȣ�m

		if (sort[1]==2||sort[1]==4)	//��1��������Ϊ��������-Һ����ʱҺ����������������
		{
			n_liquid=1;
			while (n_liquid<=n_liquid_total)
			{
				n=n_trans[2*n_liquid-1];	//Һ�������
				momentum=0;
				l_liquid=0;
				momentum=momentum+v_lright[n][2]*x_lright[n][2];	//����ͬ�����ܶȺͽ����
				l_liquid=l_liquid+x_lright[n][2];	
				n=n+1;
				while (n<=n_trans[2*n_liquid])	//��������Ҳ�����嵽�ҽ���
				{
					momentum=momentum+v_lleft[n][2]*x_lleft[n][2];	//������Ȩ�ۼ�
					l_liquid=l_liquid+x_lleft[n][2];
					n=n+1; 
				}
				n=n_trans[2*n_liquid-1];	//�����
				v_lright[n][2]=momentum/l_liquid;	//�µ��ٶ�
				n=n+1;
				while (n<=n_trans[2*n_liquid]-1)	//���������ٶȸ���ͬһҺ�������п�����
				{
					v_lright[n][2]=momentum/l_liquid;
					v_lleft[n][2]=momentum/l_liquid;
					n=n+1;
				}
				n=n_trans[2*n_liquid];	//�ҽ���
				v_lleft[n][2]=momentum/l_liquid;

				n_liquid=n_liquid+1;
			}
		}
		if (sort[1]==1||sort[1]==3||sort[1]==8)	//��1��������ΪҺ��Һ-����Һ-��-Һ����ʱ�����������
		{
			n_liquid=1;
			n=n_trans[2*n_liquid_total];	//��1��Һ���������
			momentum=0;
			l_liquid=0;
			momentum=momentum+v_lright[n][2]*x_lright[n][2];
			l_liquid=l_liquid+x_lright[n][2];
			n=n+1;
			while (n<=length)	//Һ��������ұ߿����嵽����ԭ��
			{
				momentum=momentum+v_lleft[n][2]*x_lleft[n][2];
				l_liquid=l_liquid+x_lleft[n][2];
				n=n+1;
			}
			n=1;
			while (n<=n_trans[1])	//����ԭ�㵽�ҽ���
			{
				momentum=momentum+v_lleft[n][2]*x_lleft[n][2];
				l_liquid=l_liquid+x_lleft[n][2];
				n=n+1;
			}
			n=n_trans[2*n_liquid_total];
			v_lright[n][2]=momentum/l_liquid;	//������µ��ٶ�
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

			n_liquid=2;	//�ӵ�2��Һ����ʼ��������ԭ��
			while (n_liquid<=n_liquid_total)
			{
				n=n_trans[2*n_liquid-2];	//�����
				momentum=0;
				l_liquid=0;
				momentum=momentum+v_lright[n][2]*x_lright[n][2];
				l_liquid=l_liquid+x_lright[n][2];
				n=n+1;
				while (n<=n_trans[2*n_liquid-1])	//������ұ߿����嵽�ҽ���
				{
					momentum=momentum+v_lleft[n][2]*x_lleft[n][2];
					l_liquid=l_liquid+x_lleft[n][2];
					n=n+1;
				}
				n=n_trans[2*n_liquid-2];
				v_lright[n][2]=momentum/l_liquid;	//������µ��ٶ�
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

		/*%%%%%%%%ȷ�������ݵ�״̬--������������ܶȡ��¶ȡ�ѹ��%%%%%%%%%*/
		//double c_vv[100];	//ÿ�����ݵĵ��ݱ��ȣ�J/kg-K
		//double c_pv[100];	//ÿ�����ݵĵ�ѹ���ȣ�J/kg-K
		double u_vb[100][3]; //ÿ�����ݵ����ܣ�J/kg

		if (sort[1]==2||sort[1]==4)	//��1��������Ϊ������-Һ
		{
			n_bubble=1;
			while (n_bubble<=n_bubble_total)
			{
				//c_pv[n_bubble]=-1762920+399530*T_vb[n_bubble][1]-35873*pow(T_vb[n_bubble][1],2)+1605.24*pow(T_vb[n_bubble][1],3)-35.8187*pow(T_vb[n_bubble][1],4)+0.319537*pow(T_vb[n_bubble][1],5);
				//c_vv[n_bubble]=9711.08-615.358*T_vb[n_bubble][1]+41.7627*pow(T_vb[n_bubble][1],2)-1.2972*pow(T_vb[n_bubble][1],3)+0.0168863*pow(T_vb[n_bubble][1],4);

				if (n_bubble==1)	//��һ�����ݵ������ǣ���Ϊ���������ԭ��
				{
					V_bubble[n_bubble][2]=A_v*(x_v[n_trans[2*n_bubble-1]][2]+x_v[n_trans[2*n_bubble_total]][2]+(length-1+n_trans[2*n_bubble-1]-n_trans[2*n_bubble_total])*dx);
					m_bubble[n_bubble][2]=m_bubble[n_bubble][1]+Tr_v[n_bubble]*dt+Tr_vl[n_bubble]*dt;	//��������=ԭ����+ҺĤ����+Һ�����
					ro_vb[n_bubble][2]=m_bubble[n_bubble][2]/V_bubble[n_bubble][2];	//�ܶ�=����/�����������ͨ���¶Ⱥ�ѹ��������
					/*���ݵ���������*/
					//T_vb[n_bubble][2]=1/m_bubble[n_bubble][2]*(m_bubble[n_bubble][1]*T_vb[n_bubble][1]+dt/c_vv[n_bubble]*((Tr_v[n_bubble]+Tr_vl[n_bubble])*c_pv[n_bubble]*T_vb[n_bubble][1]-P_vb[n_bubble][1]/dt*(V_bubble[n_bubble][2]-V_bubble[n_bubble][1])));	//�������ܱ仯=��������-����
					u_vb[n_bubble][1] = u_given_vt(V_bubble[n_bubble][1], m_bubble[n_bubble][1], T_vb[n_bubble][1]);//��һʱ�̵���������
					u_vb[n_bubble][2] = u_vb[n_bubble][1] + 1.0 / m_bubble[n_bubble][1] * ((Tr_v[n_bubble] + Tr_vl[n_bubble])*dt*P_vb[n_bubble][1]* V_bubble[n_bubble][1]/ m_bubble[n_bubble][1]- P_vb[n_bubble][1] * (V_bubble[n_bubble][2] - V_bubble[n_bubble][1]));//������������
					T_vb[n_bubble][2] = T_given_uv(u_vb[n_bubble][2], V_bubble[n_bubble][2], m_bubble[n_bubble][2], T_vb[n_bubble][1]);//��ʱ�̵������¶�
					P_vb[n_bubble][2]=mm*T_vb[n_bubble][2]/(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] -0.009136)-(6162* pow((1+0.1273*(1-sqrt(T_vb[n_bubble][2]/33.15))), 2))/(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] *(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] +0.009136));	//��������״̬���̸�Ϊʵ������״̬����
					P_sat= p_sat_given_T(T_vb[n_bubble][2]);	//��ȷ��
					if (P_vb[n_bubble][2]>P_sat)	//��������������Ϊ��������
					{
						P_vb[n_bubble][2]=P_sat;
						ro_vb[n_bubble][2]=ro_sat_given_T(T_vb[n_bubble][2]);	//����������״̬���̸�ΪRKS�����ʽ
						m_bubble[n_bubble][2]=ro_vb[n_bubble][2]*V_bubble[n_bubble][2];
					}
					T_sat_vb[n_bubble][2]= T_sat_given_p(P_vb[n_bubble][2]);//�������ݱ����¶�
					/*%�������ڵĸ���������г�ʼ��%*/
					n=1;	//������ԭ�㵽�ҽ���
					while (n<=n_trans[2*n_bubble-1])
					{
						P_v[n][2]=P_vb[n_bubble][2];
						T_v[n][2]=T_vb[n_bubble][2];
						T_sat[n][2]=T_sat_vb[n_bubble][2];
						ro_v[n][2]=ro_vb[n_bubble][2];
						n=n+1;
					}
					n=n_trans[2*n_bubble_total];	//������浽����ԭ��
					while (n<=length)
					{
						P_v[n][2]=P_vb[n_bubble][2];
						T_v[n][2]=T_vb[n_bubble][2];
						T_sat[n][2]=T_sat_vb[n_bubble][2];
						ro_v[n][2]=ro_vb[n_bubble][2];
						n=n+1;
					}
					
					Tr_v[n_bubble]=0;
					n=n_trans[2*n_bubble_total];	//������浽����ԭ��
					while (n<=length)
					{
						P_sat= p_sat_given_T(T_w[n][1]);
						if (P_sat>=P_v[n][1])	//�����¶ȶ�Ӧ�ı���ѹ����������ѹ��
						{
							Tr_v_sig[n]=1;	//����
						}
						else
						{
							Tr_v_sig[n]=-1;	//����
						}

						h_fg_v[n]= h_fg_given_T(T_v[n][1]);//Ϊ�������У�
						
						Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//�������������������ۼ�
						n=n+1;
					}
					n=1;
					while (n<=n_trans[2*n_bubble-1])	//������ԭ�㵽�ҽ���
					{
						P_sat= p_sat_given_T(T_w[n][1]);
						if (P_sat>=P_v[n][1])	//�����¶ȶ�Ӧ�ı���ѹ����������ѹ��
						{
							Tr_v_sig[n]=1;	//����
						}
						else
						{
							Tr_v_sig[n]=-1;
						}

						h_fg_v[n]= h_fg_given_T(T_v[n][1]);//Ϊ�������У�

						Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//�������������������ۼ�
						n=n+1;
					}
				}
				else  //�������������ԭ�������
				{
					if (n_trans[2*n_bubble-1]==n_trans[2*n_bubble-2])	//�������ҽ�����ͬһ��������
					{
						V_bubble[n_bubble][2]=A_v*x_v[n_trans[2*n_bubble-2]][2];
					}
					else  //���ҽ��治��ͬһ������
					{
						V_bubble[n_bubble][2]=A_v*(x_v[n_trans[2*n_bubble-2]][2]+x_v[n_trans[2*n_bubble-1]][2]+(n_trans[2*n_bubble-1]-n_trans[2*n_bubble-2]-1)*dx);
					}
					m_bubble[n_bubble][2]=m_bubble[n_bubble][1]+Tr_v[n_bubble]*dt+Tr_vl[n_bubble]*dt;	//�����µ�����
					ro_vb[n_bubble][2]=m_bubble[n_bubble][2]/V_bubble[n_bubble][2];	//�����µ��ܶ�
					/*���ݵ���������*/
					//T_vb[n_bubble][2]=1/m_bubble[n_bubble][2]*(m_bubble[n_bubble][1]*T_vb[n_bubble][1]+dt/c_vv[n_bubble]*((Tr_v[n_bubble]+Tr_vl[n_bubble])*c_pv[n_bubble]*T_vb[n_bubble][1]-P_vb[n_bubble][1]/dt*(V_bubble[n_bubble][2]-V_bubble[n_bubble][1])));	//�����µ��¶�
					u_vb[n_bubble][1] = u_given_vt(V_bubble[n_bubble][1], m_bubble[n_bubble][1], T_vb[n_bubble][1]);//��һʱ�̵���������
					u_vb[n_bubble][2] = u_vb[n_bubble][1] + 1.0 / m_bubble[n_bubble][1] * ((Tr_v[n_bubble] + Tr_vl[n_bubble])*dt*P_vb[n_bubble][1] * V_bubble[n_bubble][1] / m_bubble[n_bubble][1] - P_vb[n_bubble][1] * (V_bubble[n_bubble][2] - V_bubble[n_bubble][1]));//������������
					T_vb[n_bubble][2] = T_given_uv(u_vb[n_bubble][2], V_bubble[n_bubble][2], m_bubble[n_bubble][2], T_vb[n_bubble][1]);//��ʱ�̵������¶�
					P_vb[n_bubble][2]=mm*T_vb[n_bubble][2]/(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] -0.009136)-(6162* pow((1+0.1273*(1-sqrt(T_vb[n_bubble][2]/33.15))), 2))/(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] *(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] +0.009136));	//��������״̬���̸�Ϊʵ������״̬����
					P_sat= p_sat_given_T(T_vb[n_bubble][2]);	//�����¶ȶ�Ӧ�ı���ѹ��
					if (P_vb[n_bubble][2]>P_sat)	//�������ݣ�����Ϊ��������
					{
						P_vb[n_bubble][2]=P_sat;
						ro_vb[n_bubble][2]= ro_sat_given_T(T_vb[n_bubble][2]);	//����������״̬���̸�ΪRKS�����ʽ
						m_bubble[n_bubble][2]=ro_vb[n_bubble][2]*V_bubble[n_bubble][2];
					}
					T_sat_vb[n_bubble][2]= T_sat_given_p(P_vb[n_bubble][2]);//�������ݱ����¶�
					/*%�������ڵĸ���������г�ʼ��%*/
					n=n_trans[2*n_bubble-2];	//�����
					while (n<=n_trans[2*n_bubble-1])	//���ҽ���
					{
						P_v[n][2]=P_vb[n_bubble][2];
						T_v[n][2]=T_vb[n_bubble][2];
						T_sat[n][2] = T_sat_vb[n_bubble][2];
						ro_v[n][2]=ro_vb[n_bubble][2];
						n=n+1;
					}
				
					Tr_v[n_bubble]=0;
					n=n_trans[2*n_bubble-2];	//�����
					while (n<=n_trans[2*n_bubble-1])	//���ҽ���
					{
						P_sat= p_sat_given_T(T_w[n][1]);
						if (P_sat>=P_v[n][1])	//�����¶ȶ�Ӧ�ı���ѹ����������ѹ��
						{
							Tr_v_sig[n]=1;	//����
						}
						else
						{
							Tr_v_sig[n]=-1;	//����
						}
						h_fg_v[n]= h_fg_given_T(T_v[n][1]);//Ϊ�������У�
						Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//�������������������ۼ�
						n=n+1;
					}
				}
				
				/*%�����ݵ�״ֵ̬���¸���%*/
				P_vb[n_bubble][1]=P_vb[n_bubble][2];
				T_vb[n_bubble][1]=T_vb[n_bubble][2];
				T_sat_vb[n_bubble][1] = T_sat_vb[n_bubble][2];
				ro_vb[n_bubble][1]=ro_vb[n_bubble][2];
				V_bubble[n_bubble][1]=V_bubble[n_bubble][2];
				m_bubble[n_bubble][1]=m_bubble[n_bubble][2];
             	n_bubble=n_bubble+1;
			}
			/*Һ����ѹ�����Էֲ�*/
			n_bubble=1;
			while (n_bubble<=n_bubble_total)	//������ѹ������Һ��
			{
				P_l[n_trans[2*n_bubble-1]][2]=P_v[n_trans[2*n_bubble-1]][2];
				P_l[n_trans[2*n_bubble]][2]=P_v[n_trans[2*n_bubble]][2];
				n=n_trans[2*n_bubble-1]+1;	//Һ��������Ҳ������
				while (n<=n_trans[2*n_bubble]-1)	//��Һ���ҽ�����������
				{
					P_l[n][2]=P_l[n_trans[2*n_bubble]][2]-((n_trans[2*n_bubble]-n-0.5)*dx+x_lleft[n_trans[2*n_bubble]][2]/2)
						/((n_trans[2*n_bubble]-n_trans[2*n_bubble-1]-1)*dx+x_lleft[n_trans[2*n_bubble]][2]+x_lright[n_trans[2*n_bubble-1]][2])
						*(P_l[n_trans[2*n_bubble]][2]-P_l[n_trans[2*n_bubble-1]][2]);	//ÿ����������Һ���ѹ����ʾ��Һ������λ��
					n=n+1;
				}
				n_bubble=n_bubble+1;
			}
		}
		else  //��1����������Һ�忪ʼ�������������Һ��
		{
			n_bubble=1;
			while (n_bubble<=n_bubble_total)
			{
				//c_pv[n_bubble]=-1762920+399530*T_vb[n_bubble][1]-35873*pow(T_vb[n_bubble][1],2)+1605.24*pow(T_vb[n_bubble][1],3)-35.8187*pow(T_vb[n_bubble][1],4)+0.319537*pow(T_vb[n_bubble][1],5);
				//c_vv[n_bubble]=9711.08-615.358*T_vb[n_bubble][1]+41.7627*pow(T_vb[n_bubble][1],2)-1.2972*pow(T_vb[n_bubble][1],3)+0.0168863*pow(T_vb[n_bubble][1],4);
				if (n_trans[2*n_bubble]==n_trans[2*n_bubble-1])	//���������
				{
					V_bubble[n_bubble][2]=A_v*x_v[n_trans[2*n_bubble-1]][2];
				}
				else
				{
					V_bubble[n_bubble][2]=A_v*(x_v[n_trans[2*n_bubble-1]][2]+x_v[n_trans[2*n_bubble]][2]+(n_trans[2*n_bubble]-n_trans[2*n_bubble-1]-1)*dx);
				}
				m_bubble[n_bubble][2]=m_bubble[n_bubble][1]+Tr_v[n_bubble]*dt+Tr_vl[n_bubble]*dt;	//����������
				ro_vb[n_bubble][2]=m_bubble[n_bubble][2]/V_bubble[n_bubble][2];	//�������ܶ�
				/*���ݵ���������*/
				//T_vb[n_bubble][2]=1/m_bubble[n_bubble][2]*(m_bubble[n_bubble][1]*T_vb[n_bubble][1]+dt/c_vv[n_bubble]*((Tr_v[n_bubble]+Tr_vl[n_bubble])*c_pv[n_bubble]*T_vb[n_bubble][1]-P_vb[n_bubble][1]/dt*(V_bubble[n_bubble][2]-V_bubble[n_bubble][1])));	//ԭ�������̣�����
				u_vb[n_bubble][1] = u_given_vt(V_bubble[n_bubble][1], m_bubble[n_bubble][1], T_vb[n_bubble][1]);//��һʱ�̵���������
				u_vb[n_bubble][2] = u_vb[n_bubble][1] + 1.0 / m_bubble[n_bubble][1] * ((Tr_v[n_bubble] + Tr_vl[n_bubble])*dt*P_vb[n_bubble][1] * V_bubble[n_bubble][1] / m_bubble[n_bubble][1] - P_vb[n_bubble][1] * (V_bubble[n_bubble][2] - V_bubble[n_bubble][1]));//������������
				T_vb[n_bubble][2] = T_given_uv(u_vb[n_bubble][2], V_bubble[n_bubble][2], m_bubble[n_bubble][2], T_vb[n_bubble][1]);//��ʱ�̵������¶�
				P_vb[n_bubble][2]=mm*T_vb[n_bubble][2]/(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] -0.009136)-(6162* pow((1+0.1273*(1-sqrt(T_vb[n_bubble][2]/33.15))), 2))/(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] *(V_bubble[n_bubble][2] / m_bubble[n_bubble][2] +0.009136));	//��������״̬���̸�Ϊʵ������״̬����
				P_sat= p_sat_given_T(T_vb[n_bubble][2]);	//�����¶ȶ�Ӧ�ı���ѹ��
				if (P_vb[n_bubble][2]>P_sat)	//�������ݱ�����Ϊ��������
				{
					P_vb[n_bubble][2]=P_sat;
					ro_vb[n_bubble][2]= ro_sat_given_T(T_vb[n_bubble][2]);	//����������״̬���̸�ΪRKS�����ʽ
					m_bubble[n_bubble][2]=ro_vb[n_bubble][2]*V_bubble[n_bubble][2];
				}
				T_sat_vb[n_bubble][2]= T_sat_given_p(P_vb[n_bubble][2]);//�������ݱ����¶�
				/*%�������ڵĸ���������г�ʼ��%*/
				n=n_trans[2*n_bubble-1];	//���������
				while (n<=n_trans[2*n_bubble])	//���ҽ���
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
					P_sat= p_sat_given_T(T_w[n][1]);	//�����¶ȶ�Ӧ�ı���ѹ��
					if (P_sat>=P_v[n][1])	//�����¶ȶ�Ӧ�ı���ѹ����������ѹ��
					{
						Tr_v_sig[n]=1;	//����
					}
					else
					{
						Tr_v_sig[n]=-1;	//����
					}
					h_fg_v[n]= h_fg_given_T(T_v[n][1]);//Ϊ�������У�
					Tr_v[n_bubble]=Tr_v[n_bubble]+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//��������������������
					n=n+1;
				}

				/*%�����ݵ�״̬���¸���%*/
				P_vb[n_bubble][1]=P_vb[n_bubble][2];
				T_vb[n_bubble][1]=T_vb[n_bubble][2];
				T_sat_vb[n_bubble][1] = T_sat_vb[n_bubble][2];
				ro_vb[n_bubble][1]=ro_vb[n_bubble][2];
				V_bubble[n_bubble][1]=V_bubble[n_bubble][2];
				m_bubble[n_bubble][1]=m_bubble[n_bubble][2];

				n_bubble=n_bubble+1;
			}

			/*������ѹ������Һ����Һ���ڲ�ѹ�����Էֲ�*/
			n_bubble=1;
			P_l[n_trans[2*n_bubble-1]][2]=P_v[n_trans[2*n_bubble-1]][2];	//Һ���ҽ���ѹ��
			P_l[n_trans[2*n_bubble_total]][2]=P_v[n_trans[2*n_bubble_total]][2];	//Һ�������ѹ��
			n=1;
			while (n<=n_trans[2*n_bubble-1]-1)	//����ԭ�㵽�ҽ�����������
			{
				P_l[n][2]=P_l[n_trans[2*n_bubble-1]][2]-((n_trans[2*n_bubble-1]-n-0.5)*dx+x_lleft[n_trans[2*n_bubble-1]][2]/2)
					/((n_trans[2*n_bubble-1]-n_trans[2*n_bubble_total]+length-1)*dx+x_lleft[n_trans[2*n_bubble-1]][2]/2+x_lright[n_trans[2*n_bubble_total]][2]/2)
					*(P_l[n_trans[2*n_bubble-1]][2]-P_l[n_trans[2*n_bubble_total]][2]);
				n=n+1;
			}
			n=n_trans[2*n_bubble_total]+1;	//Һ�������
			while (n<=length)	//������ԭ��
			{
				P_l[n][2]=P_l[n_trans[2*n_bubble_total]][2]+((n-n_trans[2*n_bubble_total]-0.5)*dx+x_lright[n_trans[2*n_bubble_total]][2]/2)
					/((n_trans[2*n_bubble-1]-n_trans[2*n_bubble_total]+length-1)*dx+x_lleft[n_trans[2*n_bubble-1]][2]/2+x_lright[n_trans[2*n_bubble_total]][2]/2)
					*(P_l[n_trans[2*n_bubble-1]][2]-P_l[n_trans[2*n_bubble_total]][2]);
				n=n+1;
			}
			n_bubble=2;
			while (n_bubble<=n_bubble_total)
			{
				P_l[n_trans[2*n_bubble-1]][2]=P_v[n_trans[2*n_bubble-1]][2];	//Һ���ҽ���ѹ��
				P_l[n_trans[2*n_bubble-2]][2]=P_v[n_trans[2*n_bubble-2]][2];	//Һ�������ѹ��
				n=n_trans[2*n_bubble-2]+1;	//Һ��������ұ߿�����
				while (n<=n_trans[2*n_bubble-1]-1)	//��Һ���ҽ�����߿�����
				{
					P_l[n][2]=P_l[n_trans[2*n_bubble-1]][2]-((n_trans[2*n_bubble-1]-n-0.5)*dx+x_lleft[n_trans[2*n_bubble-1]][2]/2)
						/((n_trans[2*n_bubble-1]-n_trans[2*n_bubble-2]-1)*dx+x_lleft[n_trans[2*n_bubble-1]][2]/2+x_lright[n_trans[2*n_bubble-2]][2]/2)
						*(P_l[n_trans[2*n_bubble-1]][2]-P_l[n_trans[2*n_bubble-2]][2]);
					n=n+1;
				}
				n_bubble=n_bubble+1;
			}
		}

		/*****************���ݲ���******************/
		double P_vb_new;	//�²������ݵ�ѹ����Pa
		double T_vb_new;	//�²������ݵ��¶ȣ�K
		double ro_vb_new;	//�²������ݵ��ܶȣ�kg/m3
		double V_bubble_new;	//�²������ݵ������m3
		double m_bubble_new;	//�²������ݵ�������kg
		double Tr_v_new;	//�²������ݵ�����������������kg
		int n_bubble_new;	//�²������ݵ����

		/*%%��һ�����ݲ�����%%*/
		n=50;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)	//�������ɵ�1���������߸���2���������Һ��
		{
			if (generate_inter1>=generate_fre/dt)	//�������ɵ�1��ʱ�����ﵽҪ��
			{
				/*�������ɵ�1���ȶȴﵽҪ��*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)//Һ�履���¶ȴ�ȷ��������
				{
					sort[n]=6;	
					x_lleft[n][2]=(dx-l_disappear)/2;	//���ݲ������ڿ��������м�
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];	//�������¶�=�����¶�
					T_sat[n][2] = T_v[n][2];	//�����ݱ����¶�=�����¶�
					P_sat= p_sat_given_T(T_v[n][2]);	//�����ݱ���ѹ��
					P_v[n][2]=P_sat;	//������ѹ��=����ѹ��
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//����������״̬���̸�ΪRKS�����ʽ;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;	//���������ܶȺ�����������

					h_fg_v[n]= h_fg_given_T(T_v[n][1]);//Ϊ�������У�

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);	//�����¶ȶ�Ӧ�ı���ѹ��
					if (P_sat>=P_v[n][2])	//�����¶�>=�����¶ȣ������
					{
						Tr_v_sig[n]=1;	//����
					}
					else  //�����ܷ�����������
					{
						Tr_v_sig[n]=-1;	//����
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//�������������������ۼ�
					
					generate_inter1=0;	//���ݲ���֮�����ü�ʱ��
					bubble_gen_total=bubble_gen_total+1;	//�²�����������+1

					/*�����ж��µ���Һ������*/
					n=1;
					if (sort[n]==3||sort[n]==4)	//Һ-��������-Һ
					{
						transfer_point_number=1;
						n_trans[1]=1;
					}
					else if (sort[n]==8||sort[n]==6)	//Һ-��-Һ����������
					{
						transfer_point_number=2;
						n_trans[1]=1;
						n_trans[2]=1;
					}
					else  //Һ������
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
					n_trans_total=transfer_point_number;                 /*%%��¼����������%%*/ 
					n_bubble_total=n_trans_total/2;                      /*%%���ݵ�����%%*/
					n_liquid_total=n_trans_total/2;                      /*%%Һ��������%%*/

					/*��������Ž����������������ݲ�����*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//��1����������Һ�忪ʼ
					{
						n_bubble=0;
						n_bubble_new=0;
						transfer_point_number = 1;	//�ӵ�1����Һ���濪ʼ���ҵ���1��������
						while (n_bubble_new==0)	//ֱ���ҵ�������Ϊֹ
						{
							if ((transfer_point_number%2)!=0)	//����
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)	//�ÿ������������ݲ���
								{
									n_bubble_new=n_bubble;	//��1�������ݵ���ţ�����ѭ��
									sort[n_trans[transfer_point_number]]=8;	//�ÿ�������sort[n]=8
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					if (sort[1]==2||sort[1]==4)	//��1�������������忪ʼ
					{
						n_bubble=1;
						n_bubble_new=0;
						transfer_point_number = 1;
						while (n_bubble_new==0)	//ֱ���ҵ�������
						{
							if ((transfer_point_number%2)==0)	//ż��
							{
								n_bubble=n_bubble+1;
								if (sort[n_trans[transfer_point_number]]==6)	//�ÿ����������ݲ���
								{
									n_bubble_new=n_bubble;	//��1�������ݵ���ţ�����ѭ��
									sort[n_trans[transfer_point_number]]=8;	//�ÿ�������sort[n]=8
								}
							}
							transfer_point_number=transfer_point_number+1;
						}
					}

					n_bubble=n_bubble_total;
					while (n_bubble>=n_bubble_new+1)	//��1��������֮��ÿ�����ݵ��������1
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
					n_bubble=n_bubble_new;	//�����ݵĲ���
					P_vb[n_bubble][1]=P_vb_new;
					T_vb[n_bubble][1]=T_vb_new;
					T_sat_vb[n_bubble][1]= T_vb[n_bubble][1];
					ro_vb[n_bubble][1]=ro_vb_new;
					V_bubble[n_bubble][1]=V_bubble_new;
					m_bubble[n_bubble][1]=m_bubble_new;
					Tr_v[n_bubble]=Tr_v_new;
				}
			}
			else  //�������ݵ�ʱ��δ��
			{
				generate_inter1=generate_inter1+1;	//ʱ�䲽������+1
			}
		}
		else  //�������ɵ���ΧҺ�岻��
		{
			generate_inter1=generate_inter1+1;	//����������+1
		}

		/*%%�ڶ������ݲ�����%%*/
		n=104;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)	//���ݲ�������ΧҺ���㹻
		{
			if (generate_inter2>=generate_fre/dt)	//ʱ�����ﵽҪ��
			{
				/*���ȶȴﵽҪ��*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)      /*%�������ݲ���%*/
				{
					sort[n]=6;	//�������ݲ���
					x_lleft[n][2]=(dx-l_disappear)/2;
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];
					T_sat[n][2]=T_v[n][2];
					P_sat= p_sat_given_T(T_v[n][2]);
					P_v[n][2]=P_sat;
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//����������״̬���̸�ΪRKS�����ʽ;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;

					h_fg_v[n]= h_fg_given_T(T_v[n][1]);//Ϊ�������У�

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);
					if (P_sat>=P_v[n][2])	//�����ʽ�����
					{
						Tr_v_sig[n]=1;	//����
					}
					else
					{
						Tr_v_sig[n]=-1;	//����
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//��������������������
					
					generate_inter2=0;	//���ü�ʱ��
					bubble_gen_total=bubble_gen_total+1;	//�²���������������+1

					/*�����ж���Һ������λ��*/
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

					n_trans_total=transfer_point_number;	//��Һ��������
					n_bubble_total=n_trans_total/2;	//�µ���������
					n_liquid_total=n_trans_total/2;		//�µ�Һ������

					/*��������Ž����������������ݲ�����*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//��1����������Һ�忪ʼ
					{
						transfer_point_number=1;
						n_bubble=0;
						n_bubble_new=0;
						while (n_bubble_new==0)	//ֱ���ҵ�������
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

					if (sort[1]==2||sort[1]==4)	//��1�������������忪ʼ
					{
						transfer_point_number=1;
						n_bubble=1;
						n_bubble_new=0;
						while (n_bubble_new==0)	//ֱ���ҵ�������
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
					while (n_bubble>=n_bubble_new+1)	//������֮����������+1
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
			else  //ʱ��������
			{
				generate_inter2=generate_inter2+1;	//ʱ�䲽������1
			}
		}
		else  //���ݲ�������ΧҺ�岻��
		{
			generate_inter2=generate_inter2+1;	//����������+1
		}

		/*%%���������ݲ�����%%*/
		n=206;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)	//���ݲ���������Һ���㹻
		{
			if (generate_inter3>=generate_fre/dt)	//ʱ�����ﵽҪ��
			{
				/*���ȶȴﵽҪ��*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)      /*%�������ݲ���%*/
				{
					sort[n]=6;
					x_lleft[n][2]=(dx-l_disappear)/2;
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];	//�������¶�=�����¶�
					T_sat[n][2]=T_v[n][2];
					P_sat= p_sat_given_T(T_v[n][2]);
					P_v[n][2]=P_sat;
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//����������״̬���̸�ΪRKS�����ʽ;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;

					h_fg_v[n] = h_fg_given_T(T_v[n][1]);//Ϊ�������У�

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);
					if (P_sat>=P_v[n][2])	//��ʽ�����
					{
						Tr_v_sig[n]=1;	//����
					}
					else
					{
						Tr_v_sig[n]=-1;
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//��������������������
					
					generate_inter3=0;	//���ü�ʱ��
					bubble_gen_total=bubble_gen_total+1;	//�²�������������+1

					/*�����ж���Һ������λ��*/
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

					n_trans_total=transfer_point_number;                 /*%%��¼����������%%*/ 
					n_bubble_total=n_trans_total/2;                      /*%%���ݵ�����%%*/
					n_liquid_total=n_trans_total/2;                      /*%%Һ��������%%*/

					/*��������Ž����������������ݲ�����*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//��1����������Һ��ʼ
					{
						transfer_point_number=1;
						n_bubble=0;
						n_bubble_new=0;
						while (n_bubble_new==0)	//ֱ���ҵ�������
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

					if (sort[1]==2||sort[1]==4)	//��1��������������ʼ
					{
						transfer_point_number=1;
						n_bubble=1;
						n_bubble_new=0;
						while (n_bubble_new==0)	//ֱ���ҵ����������
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
					while (n_bubble>=n_bubble_new+1)	//������֮����������+1
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
			else  //ʱ��������
			{
				generate_inter3=generate_inter3+1;	//ʱ�䲽������+1
			}
		}
		else  //���ݲ���������Һ�岻��
		{
			generate_inter3=generate_inter3+1;	//����������+1
		}

		/*%%���ĸ����ݲ�����%%*/
		n=260;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)	//���ݲ���������Һ���㹻
		{
			if (generate_inter4>=generate_fre/dt)	//ʱ��������Ҫ��
			{
				/*���ȶ�����Ҫ��*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)      /*%�������ݲ���%*/
				{
					sort[n]=6;
					x_lleft[n][2]=(dx-l_disappear)/2;
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];	//�������¶�=�����¶�
					T_sat[n][2] = T_v[n][2];
					P_sat= p_sat_given_T(T_v[n][2]);
					P_v[n][2]=P_sat;
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//����������״̬���̸�ΪRKS�����ʽ;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;

					h_fg_v[n] = h_fg_given_T(T_v[n][1]);;//Ϊ�������У�

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);
					if (P_sat>=P_v[n][2])	//��ʽ�����
					{
						Tr_v_sig[n]=1;	//����
					}
					else
					{
						Tr_v_sig[n]=-1;
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//������������������
					
					generate_inter4=0;	//���ò���������
					bubble_gen_total=bubble_gen_total+1;	//�²�����������+1

					/*�����ж���Һ������λ��*/
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

					n_trans_total=transfer_point_number;                 /*%%��¼����������%%*/ 
					n_bubble_total=n_trans_total/2;                      /*%%���ݵ�����%%*/
					n_liquid_total=n_trans_total/2;                      /*%%Һ��������%%*/

					/*��������Ž����������������ݲ�����*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//��1����������Һ��ʼ
					{
						transfer_point_number=1;
						n_bubble=0;
						n_bubble_new=0;
						while (n_bubble_new==0)	//ֱ���ҵ������ݵ����
						{
							if ((transfer_point_number%2)!=0)	//����
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

					if (sort[1]==2||sort[1]==4)	//��1��������������ʼ
					{
						transfer_point_number=1;
						n_bubble=1;
						n_bubble_new=0;
						while (n_bubble_new==0)	//ֱ���ҵ�������λ��
						{
							if ((transfer_point_number%2)==0)	//ż��
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
					while (n_bubble>=n_bubble_new+1)	//������֮����������+1
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
			else  //ʱ����������Ҫ��
			{
				generate_inter4=generate_inter4+1;
			}
		}
		else  //���ݲ���������Һ�岻��
		{
			generate_inter4=generate_inter4+1;	//����������+1
		}

		/*%%��������ݲ�����%%*/
		n=470;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)
		{
			if (generate_inter5>=generate_fre/dt)
			{
				/*���ȶ�����Ҫ��*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)      /*%�������ݲ���%*/
				{
					sort[n]=6;
					x_lleft[n][2]=(dx-l_disappear)/2;
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];	//�������¶�=�����¶�
					T_sat[n][2] = T_v[n][2];
					P_sat= p_sat_given_T(T_v[n][2]);
					P_v[n][2]=P_sat;
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//����������״̬���̸�ΪRKS�����ʽ;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;

					h_fg_v[n] = h_fg_given_T(T_v[n][1]);;//Ϊ�������У�

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);
					if (P_sat>=P_v[n][2])	//��ʽ�����
					{
						Tr_v_sig[n]=1;	//����
					}
					else
					{
						Tr_v_sig[n]=-1;
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//��������������������
					
					generate_inter5=0;	//���ò���������
					bubble_gen_total=bubble_gen_total+1;	//�²�����������+1

					/*�����ж���Һ������λ��*/
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

					n_trans_total=transfer_point_number;                 /*%%��¼����������%%*/ 
					n_bubble_total=n_trans_total/2;                      /*%%���ݵ�����%%*/
					n_liquid_total=n_trans_total/2;                      /*%%Һ��������%%*/

					/*��������Ž����������������ݲ�����*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//��1����������Һ��ʼ
					{
						transfer_point_number=1;
						n_bubble=0;
						n_bubble_new=0;
						while (n_bubble_new==0)	//ֱ���ҵ����������
						{
							if ((transfer_point_number%2)!=0)	//����
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

					if (sort[1]==2||sort[1]==4)	//��1��������������ʼ
					{
						transfer_point_number=1;
						n_bubble=1;
						n_bubble_new=0;
						while (n_bubble_new==0)	//ֱ���ҵ����������
						{
							if ((transfer_point_number%2)==0)	//ż��
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
					while (n_bubble>=n_bubble_new+1)	//������֮����������+1
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
			else  //ʱ��������
			{
				generate_inter5=generate_inter5+1;	//����������+1
			}
		}
		else  //���ݲ���������Һ�岻��
		{
			generate_inter5=generate_inter5+1;	//����������+1
		}

		/*%%���������ݲ�����%%*/
		n=524;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)
		{
			if (generate_inter6>=generate_fre/dt)
			{
				/*���ȶ�����Ҫ��*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)      /*%�������ݲ���%*/
				{
					sort[n]=6;
					x_lleft[n][2]=(dx-l_disappear)/2;
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];	//�������¶�=�����¶�
					T_sat[n][2] = T_v[n][2];
					P_sat= p_sat_given_T(T_v[n][2]);
					P_v[n][2]=P_sat;
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//����������״̬���̸�ΪRKS�����ʽ;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;

					h_fg_v[n] = h_fg_given_T(T_v[n][1]);//Ϊ�������У�

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);
					if (P_sat>=P_v[n][2])	//��ʽ�����
					{
						Tr_v_sig[n]=1;	//����
					}
					else
					{
						Tr_v_sig[n]=-1;
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//��������������������

					generate_inter6=0;	//���ò���������
					bubble_gen_total=bubble_gen_total+1;	//����������+1

					/*�����ж���Һ������λ��*/
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

					n_trans_total=transfer_point_number;                 /*%%��¼����������%%*/ 
					n_bubble_total=n_trans_total/2;                      /*%%���ݵ�����%%*/
					n_liquid_total=n_trans_total/2;                      /*%%Һ��������%%*/

					/*��������Ž����������������ݲ�����*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//��1����������Һ��ʼ
					{
						transfer_point_number=1;
						n_bubble=0;
						n_bubble_new=0;
						while (n_bubble_new==0)	//ֱ���ҵ����������
						{
							if ((transfer_point_number%2)!=0)	//����
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

					if (sort[1]==2||sort[1]==4)	//��1��������������ʼ
					{
						transfer_point_number=1;
						n_bubble=1;
						n_bubble_new=0;
						while (n_bubble_new==0)	//ֱ���ҵ����������
						{
							if ((transfer_point_number%2)==0)	//ż��
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
					while (n_bubble>=n_bubble_new+1)	//������֮����������+1
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
			else  //ʱ��������
			{
				generate_inter6=generate_inter6+1;
			}
		}
		else  //���ݲ���������Һ�岻��
		{
			generate_inter6=generate_inter6+1;
		}

		/*%%���߸����ݲ�����%%*/
		n=626;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)
		{
			if (generate_inter7>=generate_fre/dt)
			{
				/*���ȶ�����Ҫ��*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)      /*%�������ݲ���%*/
				{
					sort[n]=6;
					x_lleft[n][2]=(dx-l_disappear)/2;
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];	//�������¶�=�����¶�
					T_sat[n][2] = T_v[n][2];
					P_sat= p_sat_given_T(T_v[n][2]);
					P_v[n][2]=P_sat;
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//����������״̬���̸�ΪRKS�����ʽ;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;

					h_fg_v[n] = h_fg_given_T(T_v[n][1]);//Ϊ�������У�

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);
					if (P_sat>=P_v[n][2])	//��ʽ�����
					{
						Tr_v_sig[n]=1;	//����
					}
					else
					{
						Tr_v_sig[n]=-1;
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//��������������������
					
					generate_inter7=0;	//���ò���������
					bubble_gen_total=bubble_gen_total+1;	//�²�����������+1

					/*�����ж���Һ������λ��*/
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

					n_trans_total=transfer_point_number;                 /*%%��¼����������%%*/ 
					n_bubble_total=n_trans_total/2;                      /*%%���ݵ�����%%*/
					n_liquid_total=n_trans_total/2;                      /*%%Һ��������%%*/

					/*��������������������������ݲ�����*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//��1����������Һ��ʼ
					{
						transfer_point_number=1;
						n_bubble=0;
						n_bubble_new=0;
						while (n_bubble_new==0)	//ֱ���ҵ����������
						{
							if ((transfer_point_number%2)!=0)	//����
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

					if (sort[1]==2||sort[1]==4)	//��1��������������ʼ
					{
						transfer_point_number=1;
						n_bubble=1;
						n_bubble_new=0;
						while (n_bubble_new==0)	//ֱ���ҵ����������
						{
							if ((transfer_point_number%2)==0)	//ż��
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
					while (n_bubble>=n_bubble_new+1)	//������֮����������+1
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
			else  //ʱ��������
			{
				generate_inter7=generate_inter7+1;
			}
		}
		else  //���ݲ���������Һ�岻��
		{
			generate_inter7=generate_inter7+1;
		}

		/*%%�ڰ˸����ݲ�����%%*/
		n=680;
		if (sort[n]==1&&sort[n-1]==1&&sort[n-2]==1&&sort[n+1]==1&&sort[n+2]==1)
		{
			if (generate_inter8>=generate_fre/dt)
			{
				/*���ȶ���������*/
				if (T_w[n][2]-(T_sat_given_p(P_l[n][2]))>=overheat)      /*%�������ݲ���%*/
				{
					sort[n]=6;
					x_lleft[n][2]=(dx-l_disappear)/2;
					x_lright[n][2]=(dx-l_disappear)/2;
					x_v[n][2]=l_disappear;
					T_v[n][2]=T_w[n][2];	//�������¶�=�����¶�
					T_sat[n][2] = T_v[n][2];
					P_sat= p_sat_given_T(T_v[n][2]);
					P_v[n][2]=P_sat;
					ro_v[n][2]= ro_sat_given_T(T_v[n][2]);	//����������״̬���̸�ΪRKS�����ʽ;

					P_vb_new=P_v[n][2];
					T_vb_new=T_v[n][2];
					ro_vb_new=ro_v[n][2];
					V_bubble_new=A_v*x_v[n][2];
					m_bubble_new=V_bubble_new*ro_vb_new;

					h_fg_v[n] = h_fg_given_T(T_v[n][1]);//Ϊ�������У�

					Tr_v_new=0;
					P_sat= p_sat_given_T(T_w[n][1]);
					if (P_sat>=P_v[n][2])	//��ʽ�����
					{
						Tr_v_sig[n]=1;	//����
					}
					else
					{
						Tr_v_sig[n]=-1;
					}
					Tr_v_new=Tr_v_new+Tr_v_sig[n]*pi*d*h_v[n]*x_v[n][2]*fabs(T_w[n][2]-T_sat[n][2])/h_fg_v[n];	//��������������������
					
					generate_inter8=0;	//���ò���������
					bubble_gen_total=bubble_gen_total+1;	//�²�����������+1

					/*�����ж���Һ������λ��*/
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

					n_trans_total=transfer_point_number;                 /*%%��¼����������%%*/ 
					n_bubble_total=n_trans_total/2;                      /*%%���ݵ�����%%*/
					n_liquid_total=n_trans_total/2;                      /*%%Һ��������%%*/

					/*��������Ž����������������ݲ�����*/
					if (sort[1]==1||sort[1]==3||sort[1]==8)	//��1����������Һ��ʼ
					{
						transfer_point_number=1;
						n_bubble=0;
						n_bubble_new=0;
						while (n_bubble_new==0)	//ֱ���ҵ����������
						{
							if ((transfer_point_number%2)!=0)	//����
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

					if (sort[1]==2||sort[1]==4)	//��1��������������ʼ
					{
						transfer_point_number=1;
						n_bubble=1;
						n_bubble_new=0;
						while (n_bubble_new==0)	//ֱ���ҵ���1������
						{
							if ((transfer_point_number%2)==0)	//ż��
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
					while (n_bubble>=n_bubble_new+1)	//������֮����������+1
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
			else  //ʱ��������
			{
				generate_inter8=generate_inter8+1;
			}
		}
		else  //���ݲ���������Һ�岻��
		{
			generate_inter8=generate_inter8+1;
		}

		/**********************�Ծ�ֵ���и���********************/
		n=1;
		/*��ʱ�̿���������Һ���ȣ��ٶȣ��¶ȣ�ѹ�����ܶ�*/
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
		/*��ʱ�̺��¶��йص�����*/
		for (n=1;n<=length;n=n+1)
		{
			if (sort[n]==1)	//Һ
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

			if (sort[n]==2)	//��
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

			if (sort[n]==3)	//Һ-��
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

			if (sort[n]==4)	//��-Һ
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
			if (sort[n]==8)	//Һ-��-Һ
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
		/*��ʱ�̵�Һ�嵥�໻��ϵ��*/
		for(n=1;n<=length;n=n+1)
		{
			if (sort[n]==1||sort[n]==3||sort[n]==8)	//��Һ�忪ʼ�����Һ��Ļ���ϵ��
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
			else  //�����忪ʼ
			{
				Re_lleft[n]=0;
				h_lleft[n]=0;
			}
		}
		for(n=1;n<=length;n=n+1)
		{
			if (sort[n]==1||sort[n]==4||sort[n]==8)	//�Ҳ�Һ��Ļ���ϵ��
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
			else  //�Ҳ�û��Һ��
			{
				Re_lright[n]=0;
				h_lright[n]=0;
			}
		}
		/*��ʱ�̵�Һ����ڴ���ϵ��*/
		for (n=1;n<=99;n=n+1)	//ÿ��Һ�������99������ʼ��
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
			h_b[n]=0; //���ڻ���ϵ����ʼ��
		}
		/****����Һ�����ҽ����λ��������ͷ����������������������ƶ�����****/
		/**��1��������Ϊ��������-Һ**/
		if (sort[1]==2||sort[1]==4)
		{
			n_liquid=1;
			while (n_liquid<=n_liquid_total)
			{
				/*Һ���������������1*/
				if (n_trans[2*n_liquid-1]>=n_heat_start1&&n_trans[2*n_liquid-1]<=n_heat_end1)
				{
					/*�ҽ��泬��������1*/
					if (n_trans[2*n_liquid]>n_heat_end1)
					{
						n=n_trans[2*n_liquid-1];	//��Һ������浽������1�յ�
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
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//�������ڴ��ʷ�����������
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//�����������ƶ��ľ���
						/*�ҽ�����������2*/
						if ((n_trans[2*n_liquid]>=n_heat_start2)&&(n_trans[2*n_liquid]<=n_heat_end2))
						{
							n=n_heat_start2;	//��������2��㵽�ҽ�����������
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
									Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ��ҽ���
									Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
									if (sort[n]==1)
									{
										Tr_vl_left_each[n]=Tr_vl_right_each[n];
									}
								}
								n=n+1;
							}
							n=n_trans[2*n_liquid];	//�ҽ��浥������
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
							if (n_liquid==n_liquid_total)//�ҽ��洫�ʷ�����Ҳ�����
							{
								Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];	
							}
							else
							{
								Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];
							}
							dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ�������ڴ����ƶ�����
						}
						/*�ҽ��泬��������2*/
						if (n_trans[2*n_liquid]>n_heat_end2)
						{
							n=n_heat_start2;	//������2��㵽�յ�
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
									Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ��ҽ���
									Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
									if (sort[n]==1)
									{
										Tr_vl_left_each[n]=Tr_vl_right_each[n];
									}
								}
								n=n+1;
							}
							if (n_liquid==n_liquid_total)	//�ҽ�����ڴ��ʷ�����Ҳ�����
							{
								Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];
							}
							else
							{
								Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];
							}
							dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ���������ƶ��ľ���
						}
					}
					/*�ҽ�����������1*/
					else
					{
						n=n_trans[2*n_liquid-1];	//������浽�ҽ�����߿�����
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
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						n=n_trans[2*n_liquid];	//�ҽ��浥������
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
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//�������ڴ���һ�������������
						if (n_liquid==n_liquid_total)	//�������ڴ���һ�������Ҳ�����
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lleft[n_liquid]/2;
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lleft[n_liquid]/2;
						}
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//�����������ƶ��ľ���
					}	
				}
				/*Һ���������������2*/
				else if (n_trans[2*n_liquid-1]>=n_heat_start2&&n_trans[2*n_liquid-1]<=n_heat_end2)
				{
					/*�ҽ��泬��������2*/
					if (n_trans[2*n_liquid]>n_heat_end2)
					{
						n=n_trans[2*n_liquid-1];	//������浽������2�յ�
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
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//�������ڴ��ʷ�����������
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//���������ڴ����ƶ��ľ���
					}
					/*�ҽ�����������2*/
					else
					{
						n=n_trans[2*n_liquid-1];	//������浽�ҽ�����߿�����
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
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						n=n_trans[2*n_liquid];	//�ҽ��浥������
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
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//�������ڴ��ʷ���һ����������
						if (n_liquid==n_liquid_total)	//�������ڴ��ʷ���һ����Ҳ�����
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lleft[n_liquid]/2;
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lleft[n_liquid]/2;
						}
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//���������ڴ����ƶ��ľ��룬�ҽ��治�ƶ�������
					}
				}
				/*Һ���������0��������1���֮��*/
				else if (n_trans[2*n_liquid-1]>=1&&n_trans[2*n_liquid-1]<n_heat_start1)
				{
					/*�ҽ�����������1*/
					if (n_trans[2*n_liquid]>=n_heat_start1&&n_trans[2*n_liquid]<=n_heat_end1)
					{
						n=n_heat_start1;	//��������1��㵽�ҽ�����߿�����
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
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ��ҽ���
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;	
						}
						n=n_trans[2*n_liquid];	//�ҽ��浥������
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
						if (n_liquid==n_liquid_total)	//�ҽ�����ڴ��ʷ��䵽�Ҳ�����
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];
						}
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ���������ƶ��ľ���
					}
					/*�ҽ�����������1֮��������2֮ǰ*/
					else if (n_trans[2*n_liquid]>n_heat_end1&&n_trans[2*n_liquid]<n_heat_start2)
					{
						n=n_heat_start1;	//������1��㵽�յ�
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
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//�������ڴ��ʷ���һ����������
						if (n_liquid==n_liquid_total)	//�������ڴ��ʷ���һ����Ҳ�����
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lleft[n_liquid]/2;
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lleft[n_liquid]/2;
						}
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//�����������ƶ��ľ��룬�ҽ��治�ƶ�������
					}
					/*�ҽ�����������2*/
					else if (n_trans[2*n_liquid]>=n_heat_start2&&n_trans[2*n_liquid]<=n_heat_end2)
					{
						n=n_heat_start1;//������1��㵽�յ�
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
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//�������ڴ��ʷ�����������
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//�����������ƶ��ľ���
						n=n_heat_start2;	//������2��㵽�ҽ���
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
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//���ڴ����ۼӵ��ҽ���
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						if (n_liquid==n_liquid_total)	//�ҽ�����ڴ��ʷ�����Ҳ�����
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];
						}
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ���������ƶ��ľ���
					}
					/*�ҽ��泬��������2*/
					else if (n_trans[2*n_liquid]>n_heat_end2)
					{
						n=n_heat_start1;	//������1��㵽�յ�
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
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid];	//�������ڴ��ʷ�����������
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//�����������ƶ��ľ���
						n=n_heat_start2;	//������2��㵽�յ�
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
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//���ڴ����ۼӵ��ҽ���
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						if (n_liquid==n_liquid_total)	//�ҽ�����ڴ��ʷ�����Ҳ�����
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];
						}
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ���������ƶ��ľ���
					}
				}
				/*Һ���������������1֮��������2֮ǰ*/
				else if (n_trans[2*n_liquid-1]>n_heat_end1&&n_trans[2*n_liquid-1]<n_heat_start2)
				{
					/*�ҽ�����������2*/
					if (n_trans[2*n_liquid]>=n_heat_start2&&n_trans[2*n_liquid]<=n_heat_end2)
					{
						n=n_heat_start2;	//��������2��㵽�ҽ���
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
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//���ڴ����ۼӵ��ҽ���
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						if (n_liquid==n_liquid_total)	//�ҽ�����ڴ��ʷ�����Ҳ�����
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid];
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid];
						}
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ���������ƶ��ľ���
					}
					/*�ҽ��泬��������2*/
					else if (n_trans[2*n_liquid]>n_heat_end2)
					{
						n=n_heat_start2;	//������2��㵽�յ�
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
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//���ڴ����ۼӵ��ҽ���
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid]/2;	//�ҽ�����ڴ��ʷ���һ����������
						if (n_liquid==n_liquid_total)	//�ҽ�����ڴ��ʷ���һ����Ҳ�����
						{
							Tr_vl[1]=Tr_vl[1]+Tr_vl_lright[n_liquid]/2;
						}
						else
						{
							Tr_vl[n_liquid+1]=Tr_vl[n_liquid+1]+Tr_vl_lright[n_liquid]/2;
						}
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ���������ƶ��ľ���
					}
				}
				n_liquid=n_liquid+1;
			}
		}
		
		/**��1��������ΪҺ����Һ-����Һ-��-Һ**/
		else
		{
			/**��1��Һ���������ԭ�㣬�뵥�����ǣ������ҽ���ķ��ں��ƶ��ֿ�����**/
			n_liquid=1; 
			/*��1��Һ���ҽ�����������1*/
			if (n_trans[1]>=n_heat_start1&&n_trans[1]<=n_heat_end1)
			{
				n=n_heat_start1;	//������1��㵽�ҽ���
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
						Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//���ڴ����ۼӵ��ҽ���
						Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
						if (sort[n]==1)
						{
							Tr_vl_right_each[n]=Tr_vl_left_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�ҽ�����ڴ��ʷ�����Ҳ�����
				dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ���������ƶ��ľ���
			}
			/*��1��Һ���ҽ�����������1֮��������2֮ǰ*/
			else if (n_trans[1]>n_heat_end1&&n_trans[1]<n_heat_start2)
			{
				n=n_heat_start1;	//������1��㵽�յ�
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
						Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//���ڴ����ۼӵ��ҽ���
						Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
						if (sort[n]==1)
						{
							Tr_vl_right_each[n]=Tr_vl_left_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid]/2;	//�ҽ�����ڴ��ʷ���һ����Ҳ�����
				Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lright[n_liquid]/2;	//�ҽ�����ڴ��ʷ���һ����������
				dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ���������ƶ��ľ���
			}
			/*��1��Һ���ҽ�����������2*/
			else if (n_trans[1]>=n_heat_start2&&n_trans[1]<=n_heat_end2)
			{
				n=n_heat_start1;	//������1��㵽�յ�
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
						Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//���ڴ����ۼӵ������
						Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
						if (sort[n]==1)
						{
							Tr_vl_right_each[n]=Tr_vl_left_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid];	//�������ڴ��ʷ�����������
				dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//�����������ƶ��ľ���
				n=n_heat_start2;	//������2��㵽�ҽ���
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
						Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//���ڴ����ۼӵ��ҽ���
						Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
						if (sort[n]==1)
						{
							Tr_vl_right_each[n]=Tr_vl_left_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�ҽ�����ڴ��ʷ�����Ҳ�����
				dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ���������ƶ��ľ���
			}
			
			/*��1��Һ���������������1*/
			if (n_trans[2*n_liquid_total]>=n_heat_start1&&n_trans[2*n_liquid_total]<=n_heat_end1)
			{
				n=n_trans[2*n_liquid_total];	//����浽������1�յ�
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
						Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
						Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
						if (sort[n]==1)
						{
							Tr_vl_left_each[n]=Tr_vl_right_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid];	//�������ڴ��ʷ�����������
				dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//�����������ƶ��ľ���
				n=n_heat_start2;	//������1��㵽�յ�
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
						Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ��ҽ���
						Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
						if (sort[n]==1)
						{
							Tr_vl_left_each[n]=Tr_vl_right_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�ҽ�����ڴ��ʷ�����Ҳ�����
				dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
			}
			/*��1��Һ���������������1֮��������2֮ǰ*/
			else if (n_trans[2*n_liquid_total]>n_heat_end1&&n_trans[2*n_liquid_total]<n_heat_start2)
			{
				n=n_heat_start2;	//������2��㵽�յ�
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
						Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
						Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
						if (sort[n]==1)
						{
							Tr_vl_left_each[n]=Tr_vl_right_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid]/2;	//�������ڴ��ʷ���һ����������
				Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//�������ڴ��ʷ���һ����Ҳ�����
				dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ����룬�ҽ��治�ƶ�������
			}
			/*��1��Һ���������������2*/
			else if (n_trans[2*n_liquid_total]>=n_heat_start2&&n_trans[2*n_liquid_total]<=n_heat_end2)
			{
				n=n_trans[2*n_liquid_total];	//����浽������2�յ�
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
						Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
						Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
						if (sort[n]==1)
						{
							Tr_vl_left_each[n]=Tr_vl_right_each[n];
						}
					}
					n=n+1;
				}
				Tr_vl[n_liquid_total]=Tr_vl[n_liquid_total]+Tr_vl_lleft[n_liquid];	//�������ڴ��ʷ�����������
				dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
			}
			
			/**��2�������Һ��**/
			n_liquid=2; 
			while (n_liquid<=n_liquid_total)
			{
				/*�������������1*/
				if (n_trans[2*n_liquid-2]>=n_heat_start1&&n_trans[2*n_liquid-2]<=n_heat_end1)
				{
					/*�ҽ��泬��������1*/
					if (n_trans[2*n_liquid-1]>n_heat_end1)
					{
						n=n_trans[2*n_liquid-2];	//����浽������1�յ�
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
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//�������ڴ��ʷ�����������
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
						/*�ҽ�����������2*/
						if ((n_trans[2*n_liquid-1]>=n_heat_start2)&&(n_trans[2*n_liquid-1]<=n_heat_end2))
						{
							n=n_heat_start2;	//������2��㵽�ҽ�����߿�����
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
									Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ��ҽ���
									Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
									if (sort[n]==1)
									{
										Tr_vl_left_each[n]=Tr_vl_right_each[n];
									}
								}
								n=n+1;
							}
							n=n_trans[2*n_liquid-1];	//�ҽ��浥�����ǣ�����
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
							Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�ҽ�����ڴ��ʷ��䵽�Ҳ�����
							dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
						}
						/*�ҽ��泬��������2*/
						if (n_trans[2*n_liquid-1]>n_heat_end2)
						{
							n=n_heat_start2;	//������2��㵽�յ�
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
									Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ��ҽ���
									Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
									if (sort[n]==1)
									{
										Tr_vl_left_each[n]=Tr_vl_right_each[n];
									}
								}
								n=n+1;
							}
							Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�ҽ�����ڴ��ʷ�����Ҳ�����
							dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
						}
					}
					/*�ҽ�����������1*/
					else
					{
						n=n_trans[2*n_liquid-2];	//����浽�ҽ�����߿�����
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
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						n=n_trans[2*n_liquid-1];	//�ҽ��浥������
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
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid]/2;	//�������ڴ��ʷ���һ����������
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//�������ڴ��ʷ���һ����Ҳ�����
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ����룬�ҽ��治��������
					}	
				}
				/*�������������2*/
				else if (n_trans[2*n_liquid-2]>=n_heat_start2&&n_trans[2*n_liquid-2]<=n_heat_end2)
				{
					/*�ҽ��泬��������2*/
					if (n_trans[2*n_liquid-1]>n_heat_end2)
					{
						n=n_trans[2*n_liquid-2];	//����浽������2�յ�
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
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//�������ڴ��ʷ�����������
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
					}
					/*�ҽ�����������2*/
					else
					{
						n=n_trans[2*n_liquid-2];	//����浽�ҽ�����߿�����
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
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						n=n_trans[2*n_liquid-1];	//�ҽ��浥������
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
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid]/2;	//�������ڴ��ʷ���һ����������
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//�������ڴ��ʷ���һ����Ҳ�����
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ����룬�ҽ��治�ƶ�������
					}
				}
				/*�������1��������1֮��*/
				else if (n_trans[2*n_liquid-2]>=1&&n_trans[2*n_liquid-2]<n_heat_start1)
				{
					/*�ҽ�����������1*/
					if (n_trans[2*n_liquid-1]>=n_heat_start1&&n_trans[2*n_liquid-1]<=n_heat_end1)
					{
						n=n_heat_start1;	//������1��㵽�ҽ�����߿�����
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
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ��ҽ���
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;	
						}
						n=n_trans[2*n_liquid-1];	//�ҽ��浥������
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
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�ҽ�����ڴ��ʷ�����Ҳ�����
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
					}
					/*�ҽ�����������1֮��������2֮ǰ*/
					else if (n_trans[2*n_liquid-1]>n_heat_end1&&n_trans[2*n_liquid-1]<n_heat_start2)
					{
						n=n_heat_start1;	//������1��㵽�յ�
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
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid]/2;	//�������ڴ��ʷ���һ����������
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lleft[n_liquid]/2;	//�������ڴ��ʷ���һ����Ҳ�����
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);		//������ƶ����룬�ҽ��治�ƶ�������
					}
					/*�ҽ�����������2*/
					else if (n_trans[2*n_liquid-1]>=n_heat_start2&&n_trans[2*n_liquid-1]<=n_heat_end2)
					{
						n=n_heat_start1;	//������1��㵽�յ�
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
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//�������ڴ��ʷ�����������
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
						n=n_heat_start2;	//������2��㵽�ҽ���
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
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//���ڴ����ۼӵ��ҽ���
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�ҽ�����ڴ��ʷ�����Ҳ�����
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
					}
					/*�ҽ��泬��������2*/
					else if (n_trans[2*n_liquid-1]>n_heat_end2)
					{
						n=n_heat_start1;	//������1��㵽�յ�
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
								Tr_vl_lleft[n_liquid]=Tr_vl_lleft[n_liquid]+pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];	//���ڴ����ۼӵ������
								Tr_vl_right_each[n]=pi*d*x_lright[n][1]*h_b[n]*(T_w[n][1]-T_lright[n][1])/h_fg_lright[n];
								if (sort[n]==1)
								{
									Tr_vl_left_each[n]=Tr_vl_right_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lleft[n_liquid];	//�������ڴ��ʷ�����������
						dl_vl_left[n_liquid]=dl_vl_left[n_liquid]+Tr_vl_lleft[n_liquid]*dt/(A_l*ro_l);	//������ƶ�����
						n=n_heat_start2;	//������2��㵽�յ�
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
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//���ڴ����ۼӵ��ҽ���
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�ҽ�����ڴ��ʷ�����Ҳ�����
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ��涯����
					}
				}
				/*�������������1֮��������2֮ǰ*/
				else if (n_trans[2*n_liquid-2]>n_heat_end1&&n_trans[2*n_liquid-2]<n_heat_start2)
				{
					/*�ҽ�����������2*/
					if (n_trans[2*n_liquid-1]>=n_heat_start2&&n_trans[2*n_liquid-1]<=n_heat_end2)
					{
						n=n_heat_start2;	//������2��㵽�ҽ���
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
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//���ڴ����ۼӵ��ҽ���
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid];	//�ҽ�����ڴ��ʷ�����Ҳ�����
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ�����
					}
					/*�ҽ��泬��������2*/
					else if (n_trans[2*n_liquid-1]>n_heat_end2)
					{
						n=n_heat_start2;	//������2��㵽�յ�
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
								Tr_vl_lright[n_liquid]=Tr_vl_lright[n_liquid]+pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];	//���ڴ����ۼӵ��ҽ���
								Tr_vl_left_each[n]=pi*d*x_lleft[n][1]*h_b[n]*(T_w[n][1]-T_lleft[n][1])/h_fg_lleft[n];
								if (sort[n]==1)
								{
									Tr_vl_right_each[n]=Tr_vl_left_each[n];
								}
							}
							n=n+1;
						}
						Tr_vl[n_liquid-1]=Tr_vl[n_liquid-1]+Tr_vl_lright[n_liquid]/2;	//�ҽ�����ڴ��ʷ���һ����������
						Tr_vl[n_liquid]=Tr_vl[n_liquid]+Tr_vl_lright[n_liquid]/2;	//�ҽ�����ڴ��ʷ���һ����Ҳ�����
						dl_vl_right[n_liquid]=dl_vl_right[n_liquid]+Tr_vl_lright[n_liquid]*dt/(A_l*ro_l);	//�ҽ����ƶ����룬����治�ƶ�������
					}
				}
				n_liquid=n_liquid+1;
			}
		}
		
		/*�����趨�����*/
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
		/*����Һ���ܳ���*/
		double l_liquid_total=0;	//Һ���ܳ���
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
		/*����������Һ�峤���Ƿ���֡���Ч���֡��������*/
		int special_nan=0;	//������Ч���ֻ������ı��
		for (n=1;n<=length;n=n+1)
		{
			if (_isnan(x_lleft[n][1])!=0||_finite(x_lleft[n][1])==0||_isnan(x_lright[n][1])!=0||_finite(x_lright[n][1])==0)
			{
				special_nan=1;
				break;	//��ֹѭ��for (n=1;n<=length;n=n+1)
			}
		}
		/*���������ڱ����Ƿ���֡���Ч���֡��������*/
		for (n=1;n<=length;n=n+1)
		{
			if (_isnan(T_w[n][1])!=0||_finite(T_w[n][1])==0)
			{
				special_nan=2;
				break;	//��ֹѭ��for (n=1;n<=length;n=n+1
			}
		}
		if (special_nan!=0)
		{
			fid0=fopen("state.txt","a+");	//����a+ģʽ��һ������state.txt���ļ�������state.txt���ļ�ָ���fid0
			fprintf(fid0,"����,%d,n=%d\n",special_nan,n); //��"������special_nan��n"д����fid0ָ�����ļ�
		}


		/*%%%%%%%�����%%%%%%%%*/
		if (i%print_fre==0 )	//ʱ�䲽�������������Ҫ��
		{
			kk=kk+1;   
						
			fid0=fopen("state.txt","a+");	//����a+ģʽ��һ������state.txt���ļ�������state.txt���ļ�ָ���fid0
			fprintf(fid0,"��%d�����ݼ������\n",kk);	//��"��kk�����ݼ������"д����fid0ָ�����ļ�
			fclose(fid0);	//�ر���fid0ָ�����ļ�,���ز��������0��EOF
			printf("��%d�����ݼ������\n",kk);	//�ڿ���̨���"��kk�����ݼ������"

			fid1=fopen("ro_v.txt","a+"); //����a + ģʽ��һ������ro_v.txt���ļ�������ro_v.txt���ļ�ָ���fid1
			fid2=fopen("T_w.txt","a+");	//����a + ģʽ��һ������T_w.txt���ļ�������T_w.txt���ļ�ָ���fid2
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
			
			
			fprintf(fid1,"ro_v[n]���������ǣ�%d ",kk);	//��"���������ǣ�kk"д����fid1ָ�����ļ�ro_v.txt
			fprintf(fid2,"T_w[n]���������ǣ�%d ",kk);
			fprintf(fid3,"trans_point���������ǣ�%d ",kk);
			fprintf(fid4,"P_v[n]���������ǣ�%d ",kk);
			fprintf(fid5,"Te���������ǣ�%d ",kk);
			fprintf(fid6,"x_llert[n],x_lright[n]���������ǣ�%d ",kk);
			fprintf(fid7, "Te2���������ǣ�%d ", kk);
			fprintf(fid8,"T_v[n]���������ǣ�%d ",kk);
			fprintf(fid22, "T_sat[n]���������ǣ�%d ", kk);
			fprintf(fid9,"v_lleft[n],v_lright[n]���������ǣ�%d ",kk);
			fprintf(fid10,"T_lleft[n],T_lright[n]���������ǣ�%d ",kk);
			
			fprintf(fid13,"T_vb[n_bubble]���������ǣ�%d ",kk);
			fprintf(fid23, "T_sat_vb[n_bubble]���������ǣ�%d ", kk);
			fprintf(fid14,"P_vb[n_bubble]���������ǣ�%d ",kk);
			fprintf(fid15,"ro_vb[n_bubble]���������ǣ�%d ",kk);
			fprintf(fid16,"P_l[n]���������ǣ�%d ",kk);
			fprintf(fid18,"���������ǣ�%d ",kk);
			fprintf(fid20,"h_b[n]���������ǣ�%d ",kk);
			fprintf(fid21,"h_lleft[n],h_lright[n]���������ǣ�%d ",kk);
			
			n_bubble=1;	//���ÿ�����ݵ��¶ȡ�ѹ�����ܶ�
			while (n_bubble<=n_bubble_total)
			{
				fprintf(fid13,"T_vb[%d]=%f,",n_bubble,T_vb[n_bubble][1]);	//��"T_vb[n_bubble]=T_vb[n_bubble][1],"д����fid13ָ�����ļ�T_vb.txt
				fprintf(fid23, "T_sat_vb[%d]=%f,", n_bubble, T_sat_vb[n_bubble][1]);	//��"T_sat_vb[n_bubble]=T_sat_vb[n_bubble][1],"д����fid23ָ�����ļ�T_sat_vb.txt
				fprintf(fid14,"P_vb[%d]=%f,",n_bubble,P_vb[n_bubble][1]);	//��"P_vb[n_bubble]=P_vb[n_bubble][1],"д����fid14ָ�����ļ�P_vb.txt
				fprintf(fid15,"ro_vb[%d]=%f,",n_bubble,ro_vb[n_bubble][1]);	//��"ro_vb[n_bubble]=ro_vb[n_bubble][1],"д����fid15ָ�����ļ�ro_vb.txt
				n_bubble=n_bubble+1;
			}
			fprintf(fid13,"\n");	//������д����fid13ָ�����ļ�T_vb.txt
			fprintf(fid23, "\n");	//������д����fid23ָ�����ļ�T_sat_vb.txt
			fprintf(fid14,"\n");
			fprintf(fid15,"\n");

			
			location = 1;       /*%�������꣬�൱��n*/
			while (location<=length)	//���ÿ��������Ĳ�����д��txt�ļ�
			{
				fprintf(fid1,"%f ",ro_v[location][1]);	//%f����Ŀո����ڵ���Excelʱ����
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
				
				/*д��������1��2�ı��£��Լ����ȶ��в���ѹ��*/
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

				/*д����Һ�������ڵĿ��������꣬������������Һ�峤��*/
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

			/*����n��д�뻻��*/
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
		/*д�����ȣ�Ǳ�ȣ�Һ�峤��*/
		if ((i%heat_print_fre)==0)	//ʱ�䲽�������������Ҫ��
		{
			sensible_c_total=sensible_c_total/heat_print_fre;	//ǰ���ۼӣ��˴���ʱ�䲽������ƽ��
			sensible_e_total=sensible_e_total/heat_print_fre;
			latent_c_total=latent_c_total/heat_print_fre;
			latent_e_total=latent_e_total/heat_print_fre;
			fid11=fopen("sensible.txt","a+");
			fid12=fopen("latent.txt","a+");
			fid17=fopen("liquid_total.txt","a+");
			fprintf(fid11,"ʱ��(s)��%f,sensible_c_total=%f,sensible_e_total=%f\n",i/heat_print_fre*0.001,sensible_c_total,sensible_e_total);
			fprintf(fid12,"ʱ��(s)��%f,latent_c_total=%f,latent_e_total=%f\n",i/heat_print_fre*0.001,latent_c_total,latent_e_total);
			fprintf(fid17,"ʱ��(s)��%f,liquid_total=%f\n",i/heat_print_fre*0.001,l_liquid_total);
			fclose(fid11);
			fclose(fid12);
			fclose(fid17);
			sensible_c_total=0;	//����
			sensible_e_total=0;
			latent_c_total=0;
			latent_e_total=0;
		}
		else  //ʱ�䲽�����������Ҫ��
		{
			sensible_c_total=sensible_c_total+sensible_c;	//����һʱ�̵������������ۼ�
			sensible_e_total=sensible_e_total+sensible_e;
			latent_c_total=latent_c_total+latent_c;
			latent_e_total=latent_e_total+latent_e;
		}

		n=1;
		while (n<=length)
		{
			P_l[n][2]=0;	//����ѭ��������Һ��ѹ��
			n=n+1;
		}

		i=i+1;	//��һ��ʱ��ڵ�
	}
	

	fid0=fopen("state.txt","a+");	//����a+ģʽ��һ������state.txt���ļ�������state.txt���ļ�ָ���fid0
	if (codenumber==1||codenumber==2)
		{fprintf(fid0,"Out of range, Codenumber=%d, n=%d, i=%d",codenumber,n,i);}	//��codenumber,n,iд����fid0ָ�����ļ�
	else
		{fprintf(fid0,"The end\n");}	//����The end\n��д����fid0ָ�����ļ�
	fclose(fid0);	//�ر���fid0ָ�����ļ�,���ز��������0��EOF

	printf("The end");	//����Ļ���The end
	getchar();	//�Ӽ�����ն�����һ���ַ��������ǳ������ʱ��ֹ�����ʧ�����ܹۿ����н��������һ���ȴ�
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
