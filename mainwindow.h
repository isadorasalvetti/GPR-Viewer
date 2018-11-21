#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
	~MainWindow();

private slots:
	void on_action_Quit_triggered();

	void on_checkBoxFill_toggled(bool checked);

	void on_action_Open_triggered();

    //**** Curvature Buttons
    void on_GaussianCButton_released();
    void on_MeanCButton_released();
    //**** Smoothing Buttons
    void on_IteractiveSButton_released();

private:
	Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
