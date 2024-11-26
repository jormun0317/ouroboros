#include <iostream>
#include <math.h>
#include <string>
#include <fstream>

// coded by Jormun 20240621

namespace obs
{
    float pi = 3.1415926, time = 0, dt = 0.01, end = 10, ma, re;

    template <int S, int P>
    struct cubit
    {
        cubit *idx_0, *idx_1, *idx_2, *idx_3;
        int flg_0 = 0, flg_1 = 0, flg_2 = 0, flg_3 = 0;

        float spect[3][4][S]; // r d
        float value[4][P][P], value_x[4][P][P], value_y[4][P][P];
        float base[S][P][P], base_x[S][P][P], base_y[S][P][P];
        float test[S][P][P], test_x[S][P][P], test_y[S][P][P];
        float q[4][2], n[4][2], pos[2][P][P]; // quadrant normal positon

        float legendre(float ps, int k)
        {
            if (k == -1)
                return 0.0;
            else if (k == 0)
                return 1.0;
            else
                return ((2 * k - 1) * ps * legendre(ps, k - 1) - (k - 1) * legendre(ps, k - 2)) / k;
        }

        float degendre(float ps, int k)
        {
            if (k == -1)
                return 0.0;
            else if (k == 0)
                return 0.0;
            else
                return ((2 * k - 1) * (ps * degendre(ps, k - 1) + legendre(ps, k - 1)) - (k - 1) * degendre(ps, k - 2)) / k;
        }

        void convection_x(float cx[4], int p0, int p1)
        {
            float u[4], pres, alpha, nn, cxn[4], un[4], presn, alphan;
            u[0] = value[0][p0][p1], u[1] = value[1][p0][p1], u[2] = value[2][p0][p1], u[3] = value[3][p0][p1];
            pres = 0.4 * u[3] - 0.2 * u[1] * u[1] / u[0] - 0.2 * u[2] * u[2] / u[0],
            cx[0] = u[1], cx[1] = u[1] * u[1] / u[0] + pres, cx[2] = u[1] * u[2] / u[0], cx[3] = u[1] / u[0] * (u[3] + pres);
            int flg;
            if (p0 == 0)
                un[0] = idx_0->value[0][P - 1][p1], un[1] = idx_0->value[1][P - 1][p1], un[2] = idx_0->value[2][P - 1][p1], un[3] = idx_0->value[3][P - 1][p1], nn = n[1][0], flg = flg_0;
            else if (p0 == P - 1)
                un[0] = idx_1->value[0][0][p1], un[1] = idx_1->value[1][0][p1], un[2] = idx_1->value[2][0][p1], un[3] = idx_1->value[3][0][p1], nn = n[3][0], flg = flg_1;
            else if (p1 == 0)
                un[0] = idx_2->value[0][p0][P - 1], un[1] = idx_2->value[1][p0][P - 1], un[2] = idx_2->value[2][p0][P - 1], un[3] = idx_2->value[3][p0][P - 1], nn = n[2][0], flg = flg_2;
            else if (p1 == P - 1)
                un[0] = idx_3->value[0][p0][0], un[1] = idx_3->value[1][p0][0], un[2] = idx_3->value[2][p0][0], un[3] = idx_3->value[3][p0][0], nn = n[0][0], flg = flg_3;
            else
                return;
            if (flg == 0)
            {
                presn = 0.4 * un[3] - 0.2 * un[1] * un[1] / un[0] - 0.2 * un[2] * un[2] / un[0],
                cxn[0] = un[1], cxn[1] = un[1] * un[1] / un[0] + presn, cxn[2] = un[1] * un[2] / un[0], cxn[3] = un[1] / un[0] * (un[3] + presn);
                alpha = sqrt((u[1] * u[1] + u[2] * u[2]) / (u[0] * u[0])) + sqrt(1.4 * pres / u[0]), alphan = sqrt((un[1] * un[1] + un[2] * un[2]) / (un[0] * un[0])) + sqrt(1.4 * presn / un[0]),
                alpha = std::max(alpha, alphan);
                for (int d = 0; d < 4; ++d)
                    cx[d] = 0.5 * (cx[d] + cxn[d]) - 0.5 * alpha * nn * (u[d] - un[d]);
            }
            else
                return;
        }

        void convection_y(float cy[4], int p0, int p1)
        {
            float u[4], pres, alpha, nn, cyn[4], un[4], presn, alphan;
            u[0] = value[0][p0][p1], u[1] = value[1][p0][p1], u[2] = value[2][p0][p1], u[3] = value[3][p0][p1];
            pres = 0.4 * u[3] - 0.2 * u[1] * u[1] / u[0] - 0.2 * u[2] * u[2] / u[0],
            cy[0] = u[2], cy[1] = u[1] * u[2] / u[0], cy[2] = u[2] * u[2] / u[0] + pres, cy[3] = u[2] / u[0] * (u[3] + pres);
            int flg;
            if (p0 == 0)
                un[0] = idx_0->value[0][P - 1][p1], un[1] = idx_0->value[1][P - 1][p1], un[2] = idx_0->value[2][P - 1][p1], un[3] = idx_0->value[3][P - 1][p1], nn = n[1][1], flg = flg_0;
            else if (p0 == P - 1)
                un[0] = idx_1->value[0][0][p1], un[1] = idx_1->value[1][0][p1], un[2] = idx_1->value[2][0][p1], un[3] = idx_1->value[3][0][p1], nn = n[3][1], flg = flg_1;
            else if (p1 == 0)
                un[0] = idx_2->value[0][p0][P - 1], un[1] = idx_2->value[1][p0][P - 1], un[2] = idx_2->value[2][p0][P - 1], un[3] = idx_2->value[3][p0][P - 1], nn = n[2][1], flg = flg_2;
            else if (p1 == P - 1)
                un[0] = idx_3->value[0][p0][0], un[1] = idx_3->value[1][p0][0], un[2] = idx_3->value[2][p0][0], un[3] = idx_3->value[3][p0][0], nn = n[0][1], flg = flg_3;
            else
                return;
            if (flg == 0)
            {
                presn = 0.4 * un[3] - 0.2 * un[1] * un[1] / un[0] - 0.2 * un[2] * un[2] / un[0],
                cyn[0] = un[2], cyn[1] = un[1] * un[2] / un[0], cyn[2] = un[2] * un[2] / un[0] + presn, cyn[3] = un[2] / un[0] * (un[3] + presn);
                alpha = sqrt((u[1] * u[1] + u[2] * u[2]) / (u[0] * u[0])) + sqrt(1.4 * pres / u[0]), alphan = sqrt((un[1] * un[1] + un[2] * un[2]) / (un[0] * un[0])) + sqrt(1.4 * presn / un[0]),
                alpha = std::max(alpha, alphan);
                for (int d = 0; d < 4; ++d)
                    cy[d] = 0.5 * (cy[d] + cyn[d]) - 0.5 * alpha * nn * (u[d] - un[d]);
            }
            else
                return;
        }

        void spect_to_value(int r)
        {
            for (int d = 0; d < 4; ++d)
                for (int p0 = 0; p0 < P; ++p0)
                    for (int p1 = 0; p1 < P; ++p1)
                        value[d][p0][p1] = 0, value_x[d][p0][p1] = 0, value_y[d][p0][p1] = 0;
            for (int d = 0; d < 4; ++d)
                for (int p0 = 0; p0 < P; ++p0)
                    for (int p1 = 0; p1 < P; ++p1)
                        for (int s = 0; s < S; ++s)
                            value[d][p0][p1] += spect[r][d][s] * base[s][p0][p1], value_x[d][p0][p1] += spect[r][d][s] * base_x[s][p0][p1], value_y[d][p0][p1] += spect[r][d][s] * base_y[s][p0][p1];
        }

        void value_to_spect(int r)
        {
            for (int d = 0; d < 4; ++d)
                for (int s = 0; s < S; ++s)
                    spect[r][d][s] = 0;

            float cx[4], cy[4];

            for (int p0 = 0; p0 < P; ++p0)
                for (int p1 = 0; convection_x(cx, p0, p1), p1 < P; ++p1)
                    for (int d = 0; d < 4; ++d)
                        for (int s = 0; s < S; ++s)
                            spect[r][d][s] += cx[d] * test_x[s][p0][p1];

            for (int p0 = 0; p0 < P; ++p0)
                for (int p1 = 0; convection_y(cy, p0, p1), p1 < P; ++p1)
                    for (int d = 0; d < 4; ++d)
                        for (int s = 0; s < S; ++s)
                            spect[r][d][s] += cy[d] * test_y[s][p0][p1];
        }

        void caculation()
        {
            float l[4], ps[P], ws[P]; // position weight standard
            if (P == 3)
                ps[0] = -1.0, ps[1] = 0.0, ps[2] = +1.0,
                ws[0] = 0.0, ws[1] = 2.0, ws[2] = 0.0;
            else if (P == 4)
                ps[0] = -1.0, ps[1] = -0.57735026918963, ps[2] = +0.57735026918963, ps[3] = +1.0,
                ws[0] = 0.0, ws[1] = 1.0, ws[2] = 1.0, ws[3] = 0.0;
            else if (P == 5)
                ps[0] = -1.0, ps[1] = -0.77459666924148, ps[2] = 0.0, ps[3] = +0.77459666924148, ps[4] = +1.0,
                ws[0] = 0.0, ws[1] = 0.55555555555556, ws[2] = 0.88888888888889, ws[3] = 0.55555555555556, ws[4] = 0.0;
            else if (P == 6)
                ps[0] = -1.0, ps[1] = -0.86113631159405, ps[2] = -0.33998104358486, ps[3] = +0.33998104358486, ps[4] = +0.86113631159405, ps[5] = +1.0,
                ws[0] = 0.0, ws[1] = 0.34785484513745, ws[2] = 0.65214515486255, ws[3] = 0.65214515486255, ws[4] = 0.34785484513745, ws[5] = 0.0;
            else if (P == 7)
                ps[0] = -1.0, ps[1] = -0.90617984593866, ps[2] = -0.53846931010568, ps[3] = 0.0, ps[4] = +0.53846931010568, ps[5] = +0.90617984593866, ps[6] = +1.0,
                ws[0] = 0.0, ws[1] = 0.23692688505619, ws[2] = 0.47862867049937, ws[3] = 0.56888888888889, ws[4] = 0.47862867049937, ws[5] = 0.23692688505619, ws[6] = 0.0;
            else
                return;
            n[0][0] = q[0][1] - q[1][1], n[0][1] = q[1][0] - q[0][0], n[1][0] = q[1][1] - q[2][1], n[1][1] = q[2][0] - q[1][0],
            n[2][0] = q[2][1] - q[3][1], n[2][1] = q[3][0] - q[2][0], n[3][0] = q[3][1] - q[0][1], n[3][1] = q[0][0] - q[3][0];
            for (int i = 0; i < 4; ++i)
                l[i] = sqrt(n[i][0] * n[i][0] + n[i][1] * n[i][1]), n[i][0] /= l[i], n[i][1] /= l[i];

            float x0 = (q[0][0] - q[1][0] + q[2][0] - q[3][0]) / 4, x1 = (q[0][0] - q[1][0] - q[2][0] + q[3][0]) / 4, x2 = (q[0][0] + q[1][0] - q[2][0] - q[3][0]) / 4,
                  y0 = (q[0][1] - q[1][1] + q[2][1] - q[3][1]) / 4, y1 = (q[0][1] - q[1][1] - q[2][1] + q[3][1]) / 4, y2 = (q[0][1] + q[1][1] - q[2][1] - q[3][1]) / 4,
                  a = x1 * y0 - x0 * y1, b = x0 * y2 - x2 * y0, c = x1 * y2 - x2 * y1, x_xs, x_ys, y_xs, y_ys, jaco, j2, mass;

            int K = 1, pyramid = 0;
            while ((pyramid += K) < S)
                ++K;
            --K;
            for (int k = 0, s = 0; k <= K; ++k)
                for (int s0 = 0, s1; s1 = k - s0, s0 <= k; ++s, ++s0)
                    for (int p0 = 0; p0 < P; ++p0)
                        for (int p1 = 0; p1 < P; ++p1)
                            x_xs = x0 * ps[p1] + x1, x_ys = x0 * ps[p0] + x2, y_xs = y0 * ps[p1] + y1, y_ys = y0 * ps[p0] + y2,
                            jaco = a * ps[p0] + b * ps[p1] + c, j2 = jaco * jaco, mass = (2.0 * s0 + 1) * (2.0 * s1 + 1) / 4,
                            base[s][p0][p1] = legendre(ps[p0], s0) * legendre(ps[p1], s1) * mass,
                            test[s][p0][p1] = legendre(ps[p0], s0) * legendre(ps[p1], s1) * ws[p0] * ws[p1],
                            base_x[s][p0][p1] = (degendre(ps[p0], s0) * legendre(ps[p1], s1) * y_ys + legendre(ps[p0], s0) * degendre(ps[p1], s1) * (-y_xs)) / jaco * mass,
                            base_y[s][p0][p1] = (degendre(ps[p0], s0) * legendre(ps[p1], s1) * (-x_ys) + legendre(ps[p0], s0) * degendre(ps[p1], s1) * x_xs) / jaco * mass,
                            test_x[s][p0][p1] = (degendre(ps[p0], s0) * legendre(ps[p1], s1) * jaco - legendre(ps[p0], s0) * legendre(ps[p1], s1) * a) / j2 * y_ys * ws[p0] * ws[p1] +
                                                (legendre(ps[p0], s0) * degendre(ps[p1], s1) * jaco - legendre(ps[p0], s0) * legendre(ps[p1], s1) * b) / j2 * (-y_xs) * ws[p0] * ws[p1],
                            test_y[s][p0][p1] = (degendre(ps[p0], s0) * legendre(ps[p1], s1) * jaco - legendre(ps[p0], s0) * legendre(ps[p1], s1) * a) / j2 * (-x_ys) * ws[p0] * ws[p1] +
                                                (legendre(ps[p0], s0) * degendre(ps[p1], s1) * jaco - legendre(ps[p0], s0) * legendre(ps[p1], s1) * b) / j2 * x_xs * ws[p0] * ws[p1];
            for (int k = 0, s = 0; k <= K; ++k)
                for (int s0 = 0, s1; s1 = k - s0, s0 <= k; ++s, ++s0)
                    for (int p = 0; p < P; ++p)
                        jaco = a * ps[p] + b * ps[P - 1] + c,
                        test_x[s][p][P - 1] = legendre(ps[p], s0) * legendre(ps[P - 1], s1) / jaco * l[0] * n[0][0] * ws[p] / 2,
                        test_y[s][p][P - 1] = legendre(ps[p], s0) * legendre(ps[P - 1], s1) / jaco * l[0] * n[0][1] * ws[p] / 2,
                        jaco = a * ps[0] + b * ps[p] + c,
                        test_x[s][0][p] = legendre(ps[0], s0) * legendre(ps[p], s1) / jaco * l[1] * n[1][0] * ws[p] / 2,
                        test_y[s][0][p] = legendre(ps[0], s0) * legendre(ps[p], s1) / jaco * l[1] * n[1][1] * ws[p] / 2,
                        jaco = a * ps[p] + b * ps[0] + c,
                        test_x[s][p][0] = legendre(ps[p], s0) * legendre(ps[0], s1) / jaco * l[2] * n[2][0] * ws[p] / 2,
                        test_y[s][p][0] = legendre(ps[p], s0) * legendre(ps[0], s1) / jaco * l[2] * n[2][1] * ws[p] / 2,
                        jaco = a * ps[P - 1] + b * ps[p] + c,
                        test_x[s][P - 1][p] = legendre(ps[P - 1], s0) * legendre(ps[p], s1) / jaco * l[3] * n[3][0] * ws[p] / 2,
                        test_y[s][P - 1][p] = legendre(ps[P - 1], s0) * legendre(ps[p], s1) / jaco * l[3] * n[3][1] * ws[p] / 2;

            float x, y, xs, ys, xc = 5.0, yc = 5.0, xn, yn, rn, rho, u, v, temp, pres; // eddy value
            for (int p0 = 0; p0 < P; ++p0)
                for (int p1 = 0; p1 < P; ++p1)
                    xs = ps[p0], ys = ps[p1],
                    x = (q[0][0] * (1 + xs) * (1 + ys) + q[1][0] * (1 - xs) * (1 + ys) + q[2][0] * (1 - xs) * (1 - ys) + q[3][0] * (1 + xs) * (1 - ys)) / 4.0,
                    y = (q[0][1] * (1 + xs) * (1 + ys) + q[1][1] * (1 - xs) * (1 + ys) + q[2][1] * (1 - xs) * (1 - ys) + q[3][1] * (1 + xs) * (1 - ys)) / 4.0,
                    xn = x - xc, yn = y - yc, rn = xn * xn + yn * yn, u = 1.0 + 2.5 / pi * exp(0.5 * (1.0 - rn)) * (-yn), v = 1.0 + 2.5 / pi * exp(0.5 * (1.0 - rn)) * (+xn),
                    temp = 1.0 - 10.0 / (11.2 * pi * pi) * exp(1.0 - rn), rho = pow(temp, 2.5), pres = pow(rho, 1.4), pos[0][p0][p1] = x, pos[1][p0][p1] = y,
                    value[0][p0][p1] = rho, value[1][p0][p1] = rho * u, value[2][p0][p1] = rho * v, value[3][p0][p1] = 2.5 * pres + 0.5 * rho * (u * u + v * v);

            for (int d = 0; d < 4; ++d)
                for (int s = 0; s < S; ++s)
                    spect[0][d][s] = 0;
            for (int d = 0; d < 4; ++d)
                for (int s = 0; s < S; ++s)
                    for (int p0 = 0; p0 < P; ++p0)
                        for (int p1 = 0; p1 < P; ++p1)
                            spect[0][d][s] += value[d][p0][p1] * test[s][p0][p1];
            spect_to_value(0);
        }
    };

    template <int N0, int N1, int S, int P>
    struct block
    {
        cubit<S, P> cbt[N0][N1];
        float coord[N0 + 1][N1 + 1][2];
        float layer[N0 * P][N1 * P][7];

        block()
        {
            for (int n0 = 1; n0 < N0 - 1; ++n0)
                for (int n1 = 1; n1 < N1 - 1; ++n1)
                    cbt[n0][n1].idx_0 = &cbt[n0 - 1][n1], cbt[n0][n1].idx_1 = &cbt[n0 + 1][n1], cbt[n0][n1].idx_2 = &cbt[n0][n1 - 1], cbt[n0][n1].idx_3 = &cbt[n0][n1 + 1];
            for (int n1 = 0; n1 < N1; ++n1)
                cbt[0][n1].idx_0 = &cbt[N0 - 1][n1], cbt[0][n1].idx_1 = &cbt[1][n1], cbt[N0 - 1][n1].idx_0 = &cbt[N0 - 2][n1], cbt[N0 - 1][n1].idx_1 = &cbt[0][n1];
            for (int n1 = 1; n1 < N1 - 1; ++n1)
                cbt[0][n1].idx_2 = &cbt[0][n1 - 1], cbt[0][n1].idx_3 = &cbt[0][n1 + 1], cbt[N0 - 1][n1].idx_2 = &cbt[N0 - 1][n1 - 1], cbt[N0 - 1][n1].idx_3 = &cbt[N0 - 1][n1 + 1];
            for (int n0 = 0; n0 < N0; ++n0)
                cbt[n0][0].idx_2 = &cbt[n0][N1 - 1], cbt[n0][0].idx_3 = &cbt[n0][1], cbt[n0][N1 - 1].idx_2 = &cbt[n0][N1 - 2], cbt[n0][N1 - 1].idx_3 = &cbt[n0][0];
            for (int n0 = 1; n0 < N0 - 1; ++n0)
                cbt[n0][0].idx_0 = &cbt[n0 - 1][0], cbt[n0][0].idx_1 = &cbt[n0 + 1][0], cbt[n0][N1 - 1].idx_0 = &cbt[n0 - 1][N1 - 1], cbt[n0][N1 - 1].idx_1 = &cbt[n0 + 1][N1 - 1];
        }

        void map(float a0, float b0, float a1, float b1)
        {
            float h0 = (b0 - a0) / N0, h1 = (b1 - a1) / N1;
            for (int n0 = 0; n0 <= N0; ++n0)
                for (int n1 = 0; n1 <= N1; ++n1)
                    coord[n0][n1][0] = a0 + n0 * h0, coord[n0][n1][1] = a1 + n1 * h1;
        }
        void map(std::string file)
        {
            float trash;
            std::ifstream glasses(file, std::ios::in);
            glasses >> trash >> trash;
            for (int n0 = 0; n0 <= N0; ++n0)
                for (int n1 = 0; n1 <= N1; ++n1)
                    glasses >> coord[n0][n1][0] >> coord[n0][n1][1] >> trash;
        }

        void reality() /// ma re
        {
            for (int n0 = 0; n0 < N0; ++n0)
                for (int n1 = 0; n1 < N1; cbt[n0][n1].caculation(), ++n1)
                    for (int s = 0; s < 2; ++s)
                        cbt[n0][n1].q[0][s] = coord[n0 + 1][n1 + 1][s], cbt[n0][n1].q[1][s] = coord[n0][n1 + 1][s], cbt[n0][n1].q[2][s] = coord[n0][n1][s], cbt[n0][n1].q[3][s] = coord[n0 + 1][n1][s];
        }

        void boundary(int idx, int flg)
        {
            if (idx == 0)
                for (int n1 = 0; n1 < N1; ++n1)
                    cbt[0][n1].flg_0 = flg;
            else if (idx == 1)
                for (int n1 = 0; n1 < N1; ++n1)
                    cbt[N0 - 1][n1].flg_1 = flg;
            else if (idx == 2)
                for (int n0 = 0; n0 < N0; ++n0)
                    cbt[n0][0].flg_2 = flg;
            else if (idx == 3)
                for (int n0 = 0; n0 < N0; ++n0)
                    cbt[n0][N1 - 1].flg_3 = flg;
            else
                return;
        }

        void rk3(float dt)
        {
            for (int n0 = 0; n0 < N0; ++n0)
                for (int n1 = 0; n1 < N1; ++n1)
                    cbt[n0][n1].value_to_spect(1);
            for (int n0 = 0; n0 < N0; ++n0)
                for (int n1 = 0; n1 < N1; ++n1)
                    for (int d = 0; d < 4; ++d)
                        for (int s = 0; s < S; ++s)
                            cbt[n0][n1].spect[1][d][s] = cbt[n0][n1].spect[0][d][s] + cbt[n0][n1].spect[1][d][s] * dt;
            for (int n0 = 0; n0 < N0; ++n0)
                for (int n1 = 0; n1 < N1; ++n1)
                    cbt[n0][n1].spect_to_value(1);
            for (int n0 = 0; n0 < N0; ++n0)
                for (int n1 = 0; n1 < N1; ++n1)
                    cbt[n0][n1].value_to_spect(2);
            for (int n0 = 0; n0 < N0; ++n0)
                for (int n1 = 0; n1 < N1; ++n1)
                    for (int d = 0; d < 4; ++d)
                        for (int s = 0; s < S; ++s)
                            cbt[n0][n1].spect[2][d][s] = (3.0 / 4.0) * cbt[n0][n1].spect[0][d][s] + (1.0 / 4.0) * cbt[n0][n1].spect[1][d][s] + (1.0 / 4.0) * cbt[n0][n1].spect[2][d][s] * dt;
            for (int n0 = 0; n0 < N0; ++n0)
                for (int n1 = 0; n1 < N1; ++n1)
                    cbt[n0][n1].spect_to_value(2);
            for (int n0 = 0; n0 < N0; ++n0)
                for (int n1 = 0; n1 < N1; ++n1)
                    cbt[n0][n1].value_to_spect(1);
            for (int n0 = 0; n0 < N0; ++n0)
                for (int n1 = 0; n1 < N1; ++n1)
                    for (int d = 0; d < 4; ++d)
                        for (int s = 0; s < S; ++s)
                            cbt[n0][n1].spect[0][d][s] = (1.0 / 3.0) * cbt[n0][n1].spect[0][d][s] + (2.0 / 3.0) * cbt[n0][n1].spect[2][d][s] + (2.0 / 3.0) * cbt[n0][n1].spect[1][d][s] * dt;
            for (int n0 = 0; n0 < N0; ++n0)
                for (int n1 = 0; n1 < N1; ++n1)
                    cbt[n0][n1].spect_to_value(0);
        }

        void paraview(int interval)
        {
            static int timing = -1;
            ++timing;
            if (timing % interval)
                return;
            for (int n0 = 0; n0 < N0; ++n0)
                for (int n1 = 0; n1 < N1; ++n1)
                    for (int p0 = 0; p0 < P; ++p0)
                        for (int p1 = 0; p1 < P; ++p1)
                            layer[n0 * P + p0][n1 * P + p1][0] = cbt[n0][n1].pos[0][p0][p1], layer[n0 * P + p0][n1 * P + p1][1] = cbt[n0][n1].pos[1][p0][p1], layer[n0 * P + p0][n1 * P + p1][2] = 0,
                                                        layer[n0 * P + p0][n1 * P + p1][3] = cbt[n0][n1].value[0][p0][p1], layer[n0 * P + p0][n1 * P + p1][4] = cbt[n0][n1].value[1][p0][p1],
                                                        layer[n0 * P + p0][n1 * P + p1][5] = cbt[n0][n1].value[2][p0][p1], layer[n0 * P + p0][n1 * P + p1][6] = cbt[n0][n1].value[3][p0][p1];
            std::string file, title = "x,y,z,u0,u1,u2,u3";
            if (timing > 9999)
                file = "view_9999.csv";
            else if (timing > 999)
                file = "view_" + std::to_string(timing) + ".csv";
            else if (timing > 99)
                file = "view_0" + std::to_string(timing) + ".csv";
            else if (timing > 9)
                file = "view_00" + std::to_string(timing) + ".csv";
            else
                file = "view_000" + std::to_string(timing) + ".csv";
            std::ofstream pen(file, std::ios::out | std::ios::trunc);
            pen << title << std::endl;
            for (int p1 = 0; p1 < N1 * P; ++p1)
                for (int p0 = 0; p0 < N0 * P; ++p0)
                    pen << layer[p0][p1][0] << "," << layer[p0][p1][1] << "," << layer[p0][p1][2] << ","
                        << layer[p0][p1][3] << "," << layer[p0][p1][4] << "," << layer[p0][p1][5] << "," << layer[p0][p1][6] << std::endl;
        }
    };
}

int main()
{
    obs::block<10, 10, 10, 7> *blk = new obs::block<10, 10, 10, 7>;
    blk->map("hello.dat"), blk->reality(), blk->paraview(1);
    while (obs::time < obs::end)
        obs::time += obs::dt, blk->rk3(obs::dt), blk->paraview(10);
    return 0;
}