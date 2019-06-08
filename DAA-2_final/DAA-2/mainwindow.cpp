//!  Class for the GUI
/*!
  Class that holds all the GUI elements and slots/signals.
*/
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "main.h"
#include "clickablelabel.h"

static int STATUS=0;
static int Count;
static int i=0;

static Graph initial;
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{

    ui->setupUi(this);

    connect(ui->canvas, SIGNAL(clicked()), this, SLOT(mouse_pos()));
    ui->doubleSpinBox->setRange(0,1000000);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_2_clicked()
{
    Graph& x = getGraph();
    std::vector<Point> np = x.getNormalPoints();
    if(np.size()<4){
        ui->label->setText("Too few Points please try again");
        return;
    }
    if(STATUS==0)
    {
        Graph& x = getGraph();
        initial = x;
        Count = SLS::SLS_algo(ui->doubleSpinBox->value());
        ui->label->setText("SLS Algorithm Completed! Check Iterations");
        STATUS = 1;
    }else if(STATUS==1){
        i = 0;
        Graph& x= getGraph();
        x = initial;
        Count = SLS::SLS_algo(ui->doubleSpinBox->value());
        ui->label->setText("SLS Algorithm Completed with updated Cost! Check Iterations");
        STATUS = 1;

        Paint();
    }
    return;
}

void MainWindow::on_pushButton_clicked()
{
    if(i<Count){
        if(SLS::Iteration(i)){
            Paint();
        }else {
            std::cout<<"Error in Reading Iteration"<<std::endl;
        }
        i++;
    }else{
        ui->label->setText("End of Iterations");
    }
}

void MainWindow::on_pushButton_3_clicked()
{
    i-=2;
    //    std::cout<<"i value: "<<i<< std::endl;
    if(i>=0){
        if(SLS::Iteration(i)){
            Paint();
        }else {
            std::cout<<"Error in Reading Iteration"<<std::endl;
        }
        i++;
    }else{
        ui->label->setText("First Iteration");
        i=0;
    }
}


void MainWindow::mouse_pos()
{
//    std::cout<< ui->canvas->x << " " << ui->canvas->y << std::endl;

    if(STATUS ==0){;

        Graph& g = getGraph();
        g.addNormalPoint(ui->canvas->x,ui->canvas->y);
        Paint();

    }
}

void MainWindow::Paint(){

    int w = ui->canvas->width();
    int h = ui->canvas->height();
//    std::cout<<w <<" "<<h << std::endl;
    QPixmap pixmap(w,h);
    pixmap.fill(QColor("transparent"));

    QPainter painter(&pixmap);

    QPen pointPen(Qt::blue);
    pointPen.setCapStyle(Qt::RoundCap);
    pointPen.setWidth(10);
    painter.setRenderHint(QPainter::Antialiasing,true);
    painter.setPen(pointPen);

    Graph& g = getGraph();
    std::vector<Point> np = g.getNormalPoints();
    for(auto i : np){
//      std::cout<< i.x << " " << i.y << std::endl;
      painter.drawPoint(i.x,i.y);
    }

    QPen linepen(Qt::blue);
    linepen.setCapStyle(Qt::FlatCap);
    linepen.setWidth(1);
    painter.setRenderHint(QPainter::Antialiasing,true);
    painter.setPen(linepen);
//        for (int i=0;i<(int)np.size()-1;i++) {
//            painter.drawLine(np[i].x,np[i].y, np[i+1].x,np[i+1].y);

//        }
    std::vector<Line> nl = g.getNormalLines();
    for(auto i:nl){
        painter.drawLine(i.p1.x,i.p1.y,i.p2.x,i.p2.y);
    }

    QPen bpointPen(Qt::red);
    bpointPen.setCapStyle(Qt::RoundCap);
    bpointPen.setWidth(15);
    painter.setRenderHint(QPainter::Antialiasing,true);
    painter.setPen(bpointPen);

    std::vector<Point> bp = g.getBoldPoints();
    for(auto i : bp){
//      std::cout<< i.x << " " << i.y << std::endl;
      painter.drawPoint(i.x,i.y);
    }

    QPen blinepen(Qt::red);
    blinepen.setCapStyle(Qt::FlatCap);
    blinepen.setWidth(2);
    painter.setRenderHint(QPainter::Antialiasing,true);
    painter.setPen(blinepen);

    std::vector<Line> bl = g.getBoldLines();
    for(auto i:bl){
        painter.drawLine(i.p1.x,i.p1.y,i.p2.x,i.p2.y);
    }

    ui->canvas->setPixmap(pixmap);
    QString qstr = QString::fromStdString(g.Message);
    ui->label->setText(qstr);
}




void MainWindow::on_pushButton_4_clicked()
{

    reset();
}

void MainWindow::reset(){
    STATUS = 0;
    Graph g;
    Graph& x = getGraph();
    x = g;
    i=0;
    Paint();
}
