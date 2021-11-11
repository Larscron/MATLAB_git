function [lineDataMtx] = lineData(nbuses)
%LINEDATA Summary of this function goes here
%   Detailed explanation goes here
    switch nbuses
        case 6
            %           |  From  |   To  |     R     |     X     |     B/2  |  X'mer  |
            %           |  Bus   |  Bus  |    pu     |    pu     |     pu   | TAP (a) |
            lineDataMtx=[   1       2       0.0194      0.0592      0.0264      1;
                            1       5       0.0540      0.2230      0.0246      1;
                            2       3       0.0470      0.1980      0.0219      1;
                            2       4       0.0581      0.1763      0.0170      1;
                            2       5       0.0570      0.1379      0.0173      1;
                            3       4       0.0670      0.1710      0.0064      1;
                            4       5       0.0134      0.0421      0.0000      1;
                            5       6       0.0000      0.2520      0.0000      1.0730];
        case 14
            %         |  From |  To   |   R     |   X     |     B/2  |  X'mer  |
            %         |  Bus  | Bus   |  pu     |  pu     |     pu   | TAP (a) |
            lineDataMtx=[1      2       0.01938   0.05917    0.0264         1;
                         1      5       0.05403   0.22304    0.0246         1;
                         2      3       0.04699   0.19797    0.0219         1;
                         2      4       0.05811   0.17632    0.0170         1;
                         2      5       0.05695   0.17388    0.0173         1;
                         3      4       0.06701   0.17103    0.0064         1;
                         4      5       0.01335   0.04211    0.0            1;
                         4      7       0.0       0.20912    0.0        0.978;
                         4      9       0.0       0.55618    0.0        0.969;
                         5      6       0.0       0.25202    0.0        0.932;
                         6     11       0.09498   0.19890    0.0            1;
                         6     12       0.12291   0.25581    0.0            1;
                         6     13       0.06615   0.13027    0.0            1;
                         7      8       0.0       0.17615    0.0            1;
                         7      9       0.0       0.11001    0.0            1;
                         9     10       0.03181   0.08450    0.0            1;
                         9     14       0.12711   0.27038    0.0            1;
                        10     11       0.08205   0.19207    0.0            1;
                        12     13       0.22092   0.19988    0.0            1;
                        13     14       0.17093   0.34802    0.0            1];
        case 1401
            %           |  From  |   To  |     R     |     X     |     B/2  |  X'mer  |
            %           |  Bus   |  Bus  |    pu     |    pu     |     pu   | TAP (a) |
             lineDataMtx=[  1		2		0.01938		0.05917		0.0528/2	1;
                            1		5		0.05403		0.22304		0.0492/2	1;
                            2		3		0.04699		0.19797		0.0438/2	1;
                            2		4		0.05811		0.17632		0.034/2		1;
                            2		5		0.05695		0.17388		0.0346/2	1;
                            3		4		0.06701		0.17103		0.0128/2	1;
                            4		5		0.01335		0.04211		0           1;
                            4		7		0           0.20912		0           0.978;
                            4		9		0           0.55618		0           0.969;
                            5		6		0           0.25202		0           0.932;
                            6		11		0.09498		0.1989		0           1;
                            6		12		0.12291		0.25581		0           1;
                            6		13		0.06615		0.13027		0           1;
                            7		8		0           0.17615		0           1;
                            7		9		0           0.11001		0           1;
                            9		10		0.03181		0.0845		0           1;
                            9		14		0.12711		0.27038		0           1;
                            10		11		0.08205		0.19207		0           1;
                            12		13		0.22092		0.19988		0           1;
                            13		14		0.17093		0.34802		0           1];

        otherwise
            disp ('We do not have this matrix saved')
            lineDataMtx = 0;
    end
end
