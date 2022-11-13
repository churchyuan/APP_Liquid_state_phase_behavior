classdef APP_Run_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                    matlab.ui.Figure
        GridLayout                  matlab.ui.container.GridLayout
        LeftPanel                   matlab.ui.container.Panel
        Z_rmCEditField              matlab.ui.control.NumericEditField
        Z_rmCEditFieldLabel         matlab.ui.control.Label
        Z_rmAEditField              matlab.ui.control.NumericEditField
        Z_rmAEditFieldLabel         matlab.ui.control.Label
        ell_rmBstepEditField        matlab.ui.control.NumericEditField
        ell_rmBstepEditFieldLabel   matlab.ui.control.Label
        ell_rmBminEditField         matlab.ui.control.NumericEditField
        ell_rmBminEditFieldLabel    matlab.ui.control.Label
        ell_rmBmaxEditField         matlab.ui.control.NumericEditField
        ell_rmBmaxEditFieldLabel    matlab.ui.control.Label
        epsilon_rmAEditField        matlab.ui.control.NumericEditField
        epsilon_rmAEditFieldLabel   matlab.ui.control.Label
        epsilon_rmBEditField        matlab.ui.control.NumericEditField
        epsilon_rmBEditFieldLabel   matlab.ui.control.Label
        epsilon_rmACEditField       matlab.ui.control.NumericEditField
        epsilon_rmACEditFieldLabel  matlab.ui.control.Label
        etaEditField                matlab.ui.control.NumericEditField
        fLabel                      matlab.ui.control.Label
        N_1EditField                matlab.ui.control.NumericEditField
        N_1EditFieldLabel           matlab.ui.control.Label
        N_rmTEditField              matlab.ui.control.NumericEditField
        NLabel                      matlab.ui.control.Label
        RightPanel                  matlab.ui.container.Panel
        RunButton                   matlab.ui.control.StateButton
        LiquidStateTheoryCalculatorLabel  matlab.ui.control.Label
        StatusLamp                  matlab.ui.control.Lamp
        StatusLampLabel             matlab.ui.control.Label
        ClearButton                 matlab.ui.control.StateButton
        UIAxes                      matlab.ui.control.UIAxes
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end


    methods (Access = public)

        function results = startupFcn2(app, lamda_ei, Npi, N1i, epi, epAi, epMi, lbmax, lbmin, lbstep)
            global  lamda_e N1 Np ep epA epM zp_e za2
            lamda_e =lamda_ei;
            N1 = N1i;
            Np = Npi;
            ep = epi;
            epA = epAi;
            epM = epMi;
            zp_e = get(app.Z_rmAEditField).Value;
            za2 = get(app.Z_rmCEditField).Value;
            balance2_s(lbmax,lbmin,lbstep);%balance2_s(lb_max,lb_min,deltalb)
        end
    end
    
    methods (Access = private)
        
        function results = lamp_red(app)
            app.StatusLamp.Color=[1,0,0];
        end
        
        function results = lamp_green(app)
            app.StatusLamp.Color=[0,1,0];
        end
    end


    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: RunButton
        function RunButtonValueChanged(app, event)
            lamp_red(app);
%             waitfig = figure('CloseRequestFcn', '', 'WindowStyle', 'modal');
%         uicontrol('style', 'text', 'Units', 'normal', 'Pos', [0 0.5 1 0.1], ...
%                 'Horiz', 'center', 'Fontsize', 20, 'str', 'Calculating ... Please Wait');
             pause(0.2)
%             global  lamda_e N1 Np ep epA epM
            value = app.RunButton.Value;
            lamda_e = get(app.etaEditField).Value;
            if lamda_e==0
                lamda_e=1e-5;
            end
            Np = get(app.N_rmTEditField).Value;
            N1 = get(app.N_1EditField).Value;
            ep = get(app.epsilon_rmACEditField).Value;
            epA = get(app.epsilon_rmAEditField).Value;
            epM = get(app.epsilon_rmBEditField).Value;
            lbmax = get(app.ell_rmBmaxEditField).Value;
            lbmin = get(app.ell_rmBminEditField).Value;
            lbstep = get(app.ell_rmBstepEditField).Value;
            zp_e = get(app.Z_rmAEditField).Value;
            za2 = get(app.Z_rmCEditField).Value;
            startupFcn2(app, lamda_e, Np, N1, ep, epA, epM, lbmax, lbmin, lbstep)
            %             draw(Np,N1,lamda_e,ep,epM,epA)
            if exist([ 'Nt_',num2str(Np),'N1_',num2str(N1),'eta_',num2str(lamda_e),'epAC_',num2str(ep),'_epA_',num2str(epA),'_epM_',num2str(epM),'ZA_',num2str(zp_e),'ZC_',num2str(za2),'.mat'],"file")
            eval(['load(''', 'Nt_',num2str(Np),'N1_',num2str(N1),'eta_',num2str(lamda_e),'epAC_',num2str(ep),'_epA_',num2str(epA),'_epM_',num2str(epM),'ZA_',num2str(zp_e),'ZC_',num2str(za2),'.mat''',');']);
            data=www;
            p1=data(:,1);
            p3=data(:,2);
            lb=data(:,3);
            ax0=app.UIAxes;
            hold( ax0, 'on' )
            plot(app.UIAxes,(p1(end)+p3(end))/2,lb(end),'.','MarkerSize',18,'color',[1 0.5 0]);
            plot(app.UIAxes,p1,lb);
            plot(app.UIAxes,p3,lb);
%             set(app.RunButton,'Enable','on')
%            app.RunButtonValueChanged.Value=0;
            else
               warndlg(['There is no phase sparation when lB < ',num2str(lbmax),'. If you want to draw phase diagram please increase the maximum of lB'],'Remind')
            end
            app.RunButton.Value=0;
            lamp_green(app);
%             delete(waitfig);
        end

        % Value changed function: ClearButton
        function ClearButtonValueChanged(app, event)
            value = app.ClearButton.Value;
            ax0=app.UIAxes;
            cla(ax0)
           app.ClearButton.Value=0;
%          set(app.ClearButton,'Enable','on')
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {480, 480};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {220, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 640 480];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {220, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.BackgroundColor = [0.9412 0.9412 0.9412];
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create NLabel
            app.NLabel = uilabel(app.LeftPanel);
            app.NLabel.Interpreter = 'latex';
            app.NLabel.HorizontalAlignment = 'right';
            app.NLabel.Position = [56 342 25 22];
            app.NLabel.Text = '$N_{\rm T}$';

            % Create N_rmTEditField
            app.N_rmTEditField = uieditfield(app.LeftPanel, 'numeric');
            app.N_rmTEditField.ValueDisplayFormat = '%.0f';
            app.N_rmTEditField.Position = [95 342 100 22];
            app.N_rmTEditField.Value = 100000;

            % Create N_1EditFieldLabel
            app.N_1EditFieldLabel = uilabel(app.LeftPanel);
            app.N_1EditFieldLabel.Interpreter = 'latex';
            app.N_1EditFieldLabel.HorizontalAlignment = 'right';
            app.N_1EditFieldLabel.Position = [56 295 25 22];
            app.N_1EditFieldLabel.Text = '$N_1$';

            % Create N_1EditField
            app.N_1EditField = uieditfield(app.LeftPanel, 'numeric');
            app.N_1EditField.Position = [95 295 100 22];
            app.N_1EditField.Value = 10000;

            % Create fLabel
            app.fLabel = uilabel(app.LeftPanel);
            app.fLabel.Interpreter = 'latex';
            app.fLabel.HorizontalAlignment = 'right';
            app.fLabel.Position = [56 249 25 22];
            app.fLabel.Text = '$\eta$';

            % Create etaEditField
            app.etaEditField = uieditfield(app.LeftPanel, 'numeric');
            app.etaEditField.Position = [95 249 100 22];
            app.etaEditField.Value = 0.3;

            % Create epsilon_rmACEditFieldLabel
            app.epsilon_rmACEditFieldLabel = uilabel(app.LeftPanel);
            app.epsilon_rmACEditFieldLabel.Interpreter = 'latex';
            app.epsilon_rmACEditFieldLabel.HorizontalAlignment = 'right';
            app.epsilon_rmACEditFieldLabel.Position = [53 129 28 22];
            app.epsilon_rmACEditFieldLabel.Text = '$\epsilon_{\rm AC}$';

            % Create epsilon_rmACEditField
            app.epsilon_rmACEditField = uieditfield(app.LeftPanel, 'numeric');
            app.epsilon_rmACEditField.Position = [95 129 100 22];
            app.epsilon_rmACEditField.Value = 1;

            % Create epsilon_rmBEditFieldLabel
            app.epsilon_rmBEditFieldLabel = uilabel(app.LeftPanel);
            app.epsilon_rmBEditFieldLabel.Interpreter = 'latex';
            app.epsilon_rmBEditFieldLabel.HorizontalAlignment = 'right';
            app.epsilon_rmBEditFieldLabel.Position = [59 169 25 22];
            app.epsilon_rmBEditFieldLabel.Text = '$\epsilon_{\rm B}$';

            % Create epsilon_rmBEditField
            app.epsilon_rmBEditField = uieditfield(app.LeftPanel, 'numeric');
            app.epsilon_rmBEditField.Position = [95 169 100 22];
            app.epsilon_rmBEditField.Value = 0.2;

            % Create epsilon_rmAEditFieldLabel
            app.epsilon_rmAEditFieldLabel = uilabel(app.LeftPanel);
            app.epsilon_rmAEditFieldLabel.Interpreter = 'latex';
            app.epsilon_rmAEditFieldLabel.HorizontalAlignment = 'right';
            app.epsilon_rmAEditFieldLabel.Position = [56 208 25 22];
            app.epsilon_rmAEditFieldLabel.Text = '$\epsilon_{\rm A}$';

            % Create epsilon_rmAEditField
            app.epsilon_rmAEditField = uieditfield(app.LeftPanel, 'numeric');
            app.epsilon_rmAEditField.Position = [95 208 100 22];
            app.epsilon_rmAEditField.Value = 0.2;

            % Create ell_rmBmaxEditFieldLabel
            app.ell_rmBmaxEditFieldLabel = uilabel(app.LeftPanel);
            app.ell_rmBmaxEditFieldLabel.Interpreter = 'latex';
            app.ell_rmBmaxEditFieldLabel.HorizontalAlignment = 'right';
            app.ell_rmBmaxEditFieldLabel.Position = [31 435 50 22];
            app.ell_rmBmaxEditFieldLabel.Text = '$\ell_{\rm B}$ max';

            % Create ell_rmBmaxEditField
            app.ell_rmBmaxEditField = uieditfield(app.LeftPanel, 'numeric');
            app.ell_rmBmaxEditField.Position = [95 435 100 22];
            app.ell_rmBmaxEditField.Value = 7;

            % Create ell_rmBminEditFieldLabel
            app.ell_rmBminEditFieldLabel = uilabel(app.LeftPanel);
            app.ell_rmBminEditFieldLabel.Interpreter = 'latex';
            app.ell_rmBminEditFieldLabel.HorizontalAlignment = 'right';
            app.ell_rmBminEditFieldLabel.Position = [33 411 48 22];
            app.ell_rmBminEditFieldLabel.Text = '$\ell_{\rm B}$ min';

            % Create ell_rmBminEditField
            app.ell_rmBminEditField = uieditfield(app.LeftPanel, 'numeric');
            app.ell_rmBminEditField.Position = [95 411 100 22];
            app.ell_rmBminEditField.Value = 0.0001;

            % Create ell_rmBstepEditFieldLabel
            app.ell_rmBstepEditFieldLabel = uilabel(app.LeftPanel);
            app.ell_rmBstepEditFieldLabel.Interpreter = 'latex';
            app.ell_rmBstepEditFieldLabel.HorizontalAlignment = 'right';
            app.ell_rmBstepEditFieldLabel.Position = [33 386 48 22];
            app.ell_rmBstepEditFieldLabel.Text = '$\ell_{\rm B}$ step';

            % Create ell_rmBstepEditField
            app.ell_rmBstepEditField = uieditfield(app.LeftPanel, 'numeric');
            app.ell_rmBstepEditField.Position = [95 386 100 22];
            app.ell_rmBstepEditField.Value = 0.2;

            % Create Z_rmAEditFieldLabel
            app.Z_rmAEditFieldLabel = uilabel(app.LeftPanel);
            app.Z_rmAEditFieldLabel.Interpreter = 'latex';
            app.Z_rmAEditFieldLabel.HorizontalAlignment = 'right';
            app.Z_rmAEditFieldLabel.Position = [56 81 25 22];
            app.Z_rmAEditFieldLabel.Text = '$Z_{\rm A}$';

            % Create Z_rmAEditField
            app.Z_rmAEditField = uieditfield(app.LeftPanel, 'numeric');
            app.Z_rmAEditField.Position = [95 86 100 22];
            app.Z_rmAEditField.Value = -1;

            % Create Z_rmCEditFieldLabel
            app.Z_rmCEditFieldLabel = uilabel(app.LeftPanel);
            app.Z_rmCEditFieldLabel.Interpreter = 'latex';
            app.Z_rmCEditFieldLabel.HorizontalAlignment = 'right';
            app.Z_rmCEditFieldLabel.Position = [57 48 25 22];
            app.Z_rmCEditFieldLabel.Text = '$Z_{\rm C}$';

            % Create Z_rmCEditField
            app.Z_rmCEditField = uieditfield(app.LeftPanel, 'numeric');
            app.Z_rmCEditField.Position = [96 48 100 22];
            app.Z_rmCEditField.Value = 2;

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.BackgroundColor = [1 1 1];
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create UIAxes
            app.UIAxes = uiaxes(app.RightPanel);
            title(app.UIAxes, 'Phase Diagram')
            xlabel(app.UIAxes, '\rho_{p}\sigma^3')
            ylabel(app.UIAxes, 'l_{B}/\sigma')
            zlabel(app.UIAxes, 'Z')
            app.UIAxes.TickLabelInterpreter = 'latex';
            app.UIAxes.Position = [32 94 334 293];

            % Create ClearButton
            app.ClearButton = uibutton(app.RightPanel, 'state');
            app.ClearButton.ValueChangedFcn = createCallbackFcn(app, @ClearButtonValueChanged, true);
            app.ClearButton.Text = 'Clear';
            app.ClearButton.Position = [241 36 125 33];

            % Create StatusLampLabel
            app.StatusLampLabel = uilabel(app.RightPanel);
            app.StatusLampLabel.HorizontalAlignment = 'right';
            app.StatusLampLabel.Position = [7 444 40 22];
            app.StatusLampLabel.Text = 'Status';

            % Create StatusLamp
            app.StatusLamp = uilamp(app.RightPanel);
            app.StatusLamp.Position = [62 444 10 10];

            % Create LiquidStateTheoryCalculatorLabel
            app.LiquidStateTheoryCalculatorLabel = uilabel(app.RightPanel);
            app.LiquidStateTheoryCalculatorLabel.FontName = 'Times New Roman';
            app.LiquidStateTheoryCalculatorLabel.FontSize = 18;
            app.LiquidStateTheoryCalculatorLabel.FontWeight = 'bold';
            app.LiquidStateTheoryCalculatorLabel.Position = [94 392 271 61];
            app.LiquidStateTheoryCalculatorLabel.Text = 'Liquid-State Theory Calculator';

            % Create RunButton
            app.RunButton = uibutton(app.RightPanel, 'state');
            app.RunButton.ValueChangedFcn = createCallbackFcn(app, @RunButtonValueChanged, true);
            app.RunButton.Text = 'Run';
            app.RunButton.Position = [38 36 126 33];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = APP_Run_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end