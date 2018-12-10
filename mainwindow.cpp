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
    ui->openGLWidget->mesh.IteractiveSmoothing(intSteps);
    ui->openGLWidget->update();
}

void MainWindow::on_NSteps_valueChanged(int i){
    intSteps = i;
}

void MainWindow::on_BiIteractiveSButton_released(){
    ui->openGLWidget->mesh.BiIteractiveSmoothing(intSteps);
    ui->openGLWidget->update();
}

void MainWindow::on_BiNSteps_valueChanged(int i){
    biIntSteps = i;
}

void MainWindow::on_GsmtButton_released(){
    ui->openGLWidget->mesh.GlobalSmoothing(0);
    ui->openGLWidget->update();
}

void MainWindow::on_MagDetailsButton(){

}

void MainWindow::on_DHMapButton(){

}
