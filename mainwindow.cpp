#include <QFileDialog>
#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent) :
	QMainWindow(parent),
	ui(new Ui::MainWindow)
{
	ui->setupUi(this);
}

MainWindow::~MainWindow()
{
	delete ui;
}

void MainWindow::on_action_Quit_triggered()
{
	 QApplication::quit();
}

void MainWindow::on_checkBoxFill_toggled(bool checked)
{
	 ui->openGLWidget->setPolygonMode(checked);
}

void MainWindow::on_action_Open_triggered()
{
	QString filename = QFileDialog::getOpenFileName(this, tr("Open PLY"), ".", tr("*.ply"));

	ui->openGLWidget->loadMesh(filename);
}

//Curvatures

void MainWindow::on_GaussianCButton_released(){
    ui->openGLWidget->mesh.DisplayGaussianCurvature();
    ui->openGLWidget->update();
}

void MainWindow::on_MeanCButton_released(){
    ui->openGLWidget->mesh.DisplayMeanCurvature();
    ui->openGLWidget->update();
}

//Smoothing

void MainWindow::on_IteractiveSButton_released(){
    intSteps = ui->NSteps->value();
    ui->openGLWidget->mesh.IteractiveSmoothing(intSteps);
    ui->openGLWidget->update();
}

void MainWindow::on_BiIteractiveSButton_released(){
    biIntSteps = ui->BiNSteps->value();
    ui->openGLWidget->mesh.BiIteractiveSmoothing(biIntSteps);
    ui->openGLWidget->update();
}

void MainWindow::on_GsmtButton_released(){
    ui->openGLWidget->mesh.GlobalSmoothing(0);
    ui->openGLWidget->update();
}

void MainWindow::on_MagDetailsButton_released(){
    QVector3D L = QVector3D(ui->MagL1->value(), ui->MagL3->value(), ui->MagL2->value())/10;
    ui->openGLWidget->mesh.DetMagnification(L);
    ui->openGLWidget->update();
}

void MainWindow::on_DHMapButton_released(){

}
