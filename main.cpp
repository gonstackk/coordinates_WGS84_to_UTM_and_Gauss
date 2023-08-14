#include <QCoreApplication>
#include <QFile>
#include <QDebug>
#include <QString>
#include <complex>
#include <iostream>
#include <string>
#include <QtCore/qmath.h>
QVector <double> WGS_to_GK(double lon, double lat)
{
    double maj_axis = 6378245.000; // большая полуось (в метрах)
    double min_axis = 6356863.0188; // малая полуось (в метрах)
    double dis_lon = 28;
    double dis_lat = -130;
    lon += dis_lon;
    lat += dis_lat;
    double lon_rad = lon * M_PI / 180;//долгота в радианах
    double e = qPow(1 - qPow(min_axis, 2)/qPow(maj_axis, 2), 0.5);//эксцентриситет
    double e2 = qPow(qPow(maj_axis, 2)/qPow(min_axis, 2) - 1, 0.5);//второй эксцентриситет
    double k0 = 1; // масштабный коэффициент
    double E0 = 500000; // видимо, какое-то смещение на восток (в метрах)
    int zoneNumber = lon / 6 + 1;// номер зоны
    double lambda0_deg = (zoneNumber) * 6 - 3; // долгота центрального меридиана в грудусах
    double lambda0_rad = lambda0_deg * M_PI / 180; // долгота центрального меридиана в радианах
    double I = lon_rad - lambda0_rad;

    double k = (maj_axis - min_axis) / (maj_axis + min_axis); // коэффициент для S

    double a1 = zoneNumber *  qCos(lon_rad); // a1 - a8 -коэффициенты для расчетв x,y
    double a2 = zoneNumber * qSin(lon_rad) * qCos(lon_rad) / 2;
    double a3 = zoneNumber * qPow(qCos(lon_rad), 3) * (1 - qPow(qTan(lon_rad), 2) + qPow(e2, 2)) / 6;
    double a4 = zoneNumber * qSin(lon_rad) * qPow(qCos(lon_rad), 3) * (5 - qPow(qTan(lon_rad), 2) + 9 * qPow(e2, 2) + 4 * qPow(e2, 4)) / 24;
    double a5 = zoneNumber * qPow(qCos(lon_rad), 5) * (5 - 18 * qPow(qTan(lon_rad), 2) + qPow(qTan(lon_rad), 4) + 14 * qPow(e2, 2) - 58 * qPow(e2, 2)) * qPow(qTan(lon_rad), 2) / 120;
    double a6 = zoneNumber * qSin(lon_rad) * qPow(qCos(lon_rad), 5) * (61 - 58 * qPow(qTan(lon_rad), 2) + qPow(qTan(lon_rad), 4) + 270 * qPow(e2, 2) - 330 * qPow(e2, 2)) * qPow(qTan(lon_rad), 2) / 720;
    double a7 = zoneNumber * qPow(qCos(lon_rad), 7) * (61 - 479 * qPow(qTan(lon_rad), 2) + 179 * qPow(qTan(lon_rad), 4) - qPow(qTan(lon_rad), 6)) / 5040;
    double a8 = zoneNumber * qSin(lon_rad) * qPow(qCos(lon_rad), 7) * (1385 - 3111 * qPow(qTan(lon_rad), 2) + 543 * qPow(qTan(lon_rad), 4) - qPow(qTan(lon_rad), 6)) / 40320;

    double S = (maj_axis / (1 + k)) * ((1 + qPow(k, 2) / 4 + qPow(k, 4) / 64) * lon_rad - (3 * k / 2 - 3 * qPow(k, 3) / 16) * qSin(2 * lon_rad) +
                                       (15 * qPow(k, 2) / 16 - 15 * qPow(k, 4) / 16) * qSin(4 * lon_rad) - 35 * k * qSin(6 * lon_rad) / 48);// Начальное значение абсциссы (метры)

    //double N0 = 0;

    /*if (lat < 0) // к чему это условие? Мб север-юг? похоже на то
    {
        N0 = 10000000; // сдвиг по северу (в метрах)
    }*/
    double x = S + a2 * qPow(I, 2) + a4 * qPow(I, 4) + a6 * qPow(I, 6) + a8 * qPow(I, 8);
    double  y = a1 * I + a3 * qPow(I, 3) + a5 * qPow(I, 5) + a7 * qPow(I, 7);;



    QVector <double> vec_coor;
    vec_coor.push_back(x);
    vec_coor.push_back(y);

    return vec_coor;

}

QVector <double> WGS_to_UTM(double lon, double lat)
{
    double maj_axis = 6378137.0000; // большая полуось (в метрах)
    double min_axis = 6356752.3142; // малая полуось (в метрах)
    double e = qPow(1 - qPow(min_axis, 2)/qPow(maj_axis, 2), 0.5);//эксцентриситет
    double k0 = 0.9996; // масштабный коэффициент
    double E0 = 500000; // видимо, какое-то смещение на восток (в метрах)
    int zoneNumber = lon/6 + 31;// номер зоны
    double lambda0_deg = (zoneNumber - 1) * 6 - 180 + 3; // долгота центрального меридиана в грудусах
    double lambda0_rad = lambda0_deg * M_PI / 180; // долгота центрального меридиана в радианах
    double N0 = 0;

    if (lat < 0) // к чему это условие? Мб север-юг? похоже на то
    {
        N0 = 10000000; // сдвиг по северу (в метрах)
    }

    double lambda = lon * M_PI / 180;// долгота в радианах
    double phi = lat * M_PI / 180;// широта в радианах
    double v = 1 / sqrt(1 - qPow(e * qSin(phi), 2));// кривизна меридиана в плоскости UTM
    double A = (lambda - lambda0_rad) * qCos(phi);
    double T = qPow(qTan(phi), 2);
    double C = qPow(e, 2) / (1 - qPow(e, 2)) * qPow(qCos(phi), 2);
    double s = (1 - qPow(e, 2) / 4 - 3 * qPow(e, 4) / 64 - 5 * qPow(e, 6)/256) * phi - (3 * qPow(e, 2) / 8 + 3 * qPow(e, 4) /32
    + 45 * qPow(e, 6) / 1024) * qSin(2 * phi) + (15 * qPow(e, 4) / 256 + 45 * qPow(e, 6) / 1024) * qSin(4 * phi) - 35 * qPow(e, 6) / 3072 * qSin(6 * phi); // арктангенс главной вертикальной кривизны

    double x = E0 + k0 * maj_axis * v * (A + (1 - T + C) * qPow(A, 3)/6 + (5 - 18 * T + T * T) * qPow(A, 5) / 120);
    double y = N0 + k0 * maj_axis * (s + v * tan(phi) * (pow(A, 2)/2 + (5 - T + 9 * C + 4 * C * C) * qPow(A, 4) / 24 + (61 - 58 * T + T * T) * qPow(A, 6) / 720));

    QVector <double> vec_coor;
    vec_coor.push_back(x);
    vec_coor.push_back(y);

    return vec_coor;

}

int main()
{
    //double longitude = 157.210009;// долгота
    //double latitude = 51.712013;// широта
    double longitude = -156.32;// долгота
    double latitude = 49.66;// широта
    QVector <double> vec_x_y_WGS;
    QVector <double> vec_x_y_GK;

    QVector <double> vecIn;
    QString str;
    int N = 0;
    int flag = 0;

    QFile inputFile("/home/anastasia/C++/test/Coordinaty/coor_2.txt");
    if (inputFile.open(QIODevice::ReadOnly))
    {
        str = inputFile.readLine();
        while (!inputFile.atEnd())
        {
           bool ok;
           QString line = inputFile.readLine();
           QStringList row = line.split("\t");

           QVector <double> values(row.size());

           for(int i = 0; i< row.size(); i++)
           {
                    values[i] = row[i].toDouble(&ok);
                    if (ok == false)
                    {
                        std::cout << "исправьте ошибку в сроке №" << N / 2 + 2 << std::endl;
                        flag = 1;
                    }
                    vecIn.push_back(values[i]);
                    N++;
           }

        }
        inputFile.close();
    }
    if (flag == 0)
    {
        for (int i =0; i < N; i += 2 )
        {
            vec_x_y_WGS = WGS_to_UTM(vecIn[i], vecIn[i + 1]);
            std::cout << "x = " << vec_x_y_WGS [0] <<  " м; " << "y = " << vec_x_y_WGS [1] <<  " м;" << std::endl;
        }

        //vec_x_y_WGS = WGS_to_UTM(longitude, latitude);
        //std::cout << "x =  " << vec_x_y_WGS [0] <<  " м; " << "y =  " << vec_x_y_WGS [1] <<  " м; " << std::endl;

        //vec_x_y_GK = WGS_to_GK(longitude,latitude);
        //std::cout << "x = " << vec_x_y_GK [0] <<  " м; " << "y = " << vec_x_y_GK [1] <<  " м; " << std::endl;
    }
    return 0;
}
