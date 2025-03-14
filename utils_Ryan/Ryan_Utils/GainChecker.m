classdef GainChecker < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                     matlab.ui.Figure
        ImpedanceohmsEditFieldLabel  matlab.ui.control.Label
        InputLabel                   matlab.ui.control.Label
        InputsLabel                  matlab.ui.control.Label
        dBmEditFieldLabel            matlab.ui.control.Label
        OutputsLabel                 matlab.ui.control.Label
        VppEditFieldLabel            matlab.ui.control.Label
        ImpedanceohmsEditField       matlab.ui.control.NumericEditField
        MaximumVppEditField          matlab.ui.control.NumericEditField
        MaximumVppEditFieldLabel     matlab.ui.control.Label
        AttenuationdBEditField       matlab.ui.control.NumericEditField
        AttenuationdBEditFieldLabel  matlab.ui.control.Label
        GaindBEditField              matlab.ui.control.NumericEditField
        GaindBEditFieldLabel         matlab.ui.control.Label
        dBmEditField                 matlab.ui.control.NumericEditField
        VppEditField                 matlab.ui.control.NumericEditField
        TextArea_2                   matlab.ui.control.TextArea
        TextArea                     matlab.ui.control.TextArea
        VppTextArea                  matlab.ui.control.TextArea
        VppTextAreaLabel             matlab.ui.control.Label
        VrmsTextArea                 matlab.ui.control.TextArea
        VrmsTextAreaLabel            matlab.ui.control.Label
        PmWTextArea                  matlab.ui.control.TextArea
        PmWTextAreaLabel             matlab.ui.control.Label
        PdBmTextArea                 matlab.ui.control.TextArea
        OutputV_ppLabel              matlab.ui.control.Label
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: AttenuationdBEditField, GaindBEditField, 
        % ...and 3 other components
        function GainEditFieldValueChanged(app, event)

            R = app.ImpedanceohmsEditField.Value;
            gain = app.GaindBEditField.Value;
            attenuation = app.AttenuationdBEditField.Value;
            tag = event.Source.Tag;

        switch tag
            case "Vpp"
                inputVpp = app.VppEditField.Value;
                inputVrms = inputVpp/(2 * sqrt(2));
                inputPowermw = 1e3 * inputVrms^2 ./ R;
                inputPowerdbm = 10 * log10(inputPowermw);
                app.dBmEditField.Value = inputPowerdbm;

            case "dBm"
                inputPowerdbm = app.dBmEditField.Value;
                inputPowermw = 10^(inputPowerdbm / 10);
                inputVrms = sqrt(inputPowermw * R /1e3);
                inputVpp = 2 * sqrt(2) * inputVrms;
                app.VppEditField.Value = inputVpp;

            otherwise
                inputVpp = app.VppEditField.Value;
                inputVrms = inputVpp/(2 * sqrt(2));
                inputPowermw = 1e3 * inputVrms^2 ./ R;
                inputPowerdbm = 10 * log10(inputPowermw);
                app.dBmEditField.Value = inputPowerdbm;


        end



           
            outputPowerdbm = inputPowerdbm + gain - attenuation;
            outputPowermw = 10 ^ (outputPowerdbm / 10);
            outputVrms = sqrt(outputPowermw * R / 1e3);
            outputVpp = 2 * sqrt(2) * outputVrms;
            

            
            app.PdBmTextArea.Value = num2str(outputPowerdbm, '%.1f');
            app.PmWTextArea.Value = num2str(outputPowermw, '%.3f');
            app.VrmsTextArea.Value = num2str(outputVrms, '%.3f');
            app.VppTextArea.Value = num2str(outputVpp, '%.3f');
            
            
            
            maxVpp = app.MaximumVppEditField.Value;

            minAttenuation = (10*log((125*inputVpp^2)/R) - 10*log((125*maxVpp^2)/R) + gain*log(10))/log(10);

            
            maxAttenuationRemove = attenuation - minAttenuation;

            app.TextArea.Value = "To stay under " + maxVpp + " Volts, you need a minimum of " + num2str(minAttenuation, '%.1f') + " dB of attenuation." ;
            
            if maxAttenuationRemove < 0
                app.TextArea_2.Value = "Voltage exceeds maximum value. ";

            else
                
                app.TextArea_2.Value = "You can remove a MAXIMUM of " + num2str(maxAttenuationRemove, '%.1f') + " dB of Attenuation";
            
            end
            
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 687 455];
            app.UIFigure.Name = 'MATLAB App';

            % Create OutputV_ppLabel
            app.OutputV_ppLabel = uilabel(app.UIFigure);
            app.OutputV_ppLabel.HorizontalAlignment = 'right';
            app.OutputV_ppLabel.Interpreter = 'latex';
            app.OutputV_ppLabel.Position = [30 188 78 23];
            app.OutputV_ppLabel.Text = 'P (dBm)';

            % Create PdBmTextArea
            app.PdBmTextArea = uitextarea(app.UIFigure);
            app.PdBmTextArea.Editable = 'off';
            app.PdBmTextArea.Placeholder = '0';
            app.PdBmTextArea.Position = [123 188 91 25];

            % Create PmWTextAreaLabel
            app.PmWTextAreaLabel = uilabel(app.UIFigure);
            app.PmWTextAreaLabel.HorizontalAlignment = 'right';
            app.PmWTextAreaLabel.Interpreter = 'latex';
            app.PmWTextAreaLabel.Position = [29 151 78 23];
            app.PmWTextAreaLabel.Text = 'P (mW)';

            % Create PmWTextArea
            app.PmWTextArea = uitextarea(app.UIFigure);
            app.PmWTextArea.Editable = 'off';
            app.PmWTextArea.Placeholder = '0';
            app.PmWTextArea.Position = [122 151 91 25];

            % Create VrmsTextAreaLabel
            app.VrmsTextAreaLabel = uilabel(app.UIFigure);
            app.VrmsTextAreaLabel.HorizontalAlignment = 'right';
            app.VrmsTextAreaLabel.Interpreter = 'latex';
            app.VrmsTextAreaLabel.Position = [29 111 78 23];
            app.VrmsTextAreaLabel.Text = 'Vrms';

            % Create VrmsTextArea
            app.VrmsTextArea = uitextarea(app.UIFigure);
            app.VrmsTextArea.Editable = 'off';
            app.VrmsTextArea.Placeholder = '0';
            app.VrmsTextArea.Position = [122 111 91 25];

            % Create VppTextAreaLabel
            app.VppTextAreaLabel = uilabel(app.UIFigure);
            app.VppTextAreaLabel.HorizontalAlignment = 'right';
            app.VppTextAreaLabel.Interpreter = 'latex';
            app.VppTextAreaLabel.Position = [29 70 78 23];
            app.VppTextAreaLabel.Text = 'Vpp';

            % Create VppTextArea
            app.VppTextArea = uitextarea(app.UIFigure);
            app.VppTextArea.Editable = 'off';
            app.VppTextArea.Placeholder = '0';
            app.VppTextArea.Position = [122 70 91 25];

            % Create TextArea
            app.TextArea = uitextarea(app.UIFigure);
            app.TextArea.Editable = 'off';
            app.TextArea.Position = [392 153 198 60];

            % Create TextArea_2
            app.TextArea_2 = uitextarea(app.UIFigure);
            app.TextArea_2.Editable = 'off';
            app.TextArea_2.Position = [394 76 196 60];

            % Create VppEditField
            app.VppEditField = uieditfield(app.UIFigure, 'numeric');
            app.VppEditField.ValueChangedFcn = createCallbackFcn(app, @GainEditFieldValueChanged, true);
            app.VppEditField.Tag = 'Vpp';
            app.VppEditField.Placeholder = '0';
            app.VppEditField.Position = [114 374 44 22];

            % Create dBmEditField
            app.dBmEditField = uieditfield(app.UIFigure, 'numeric');
            app.dBmEditField.ValueChangedFcn = createCallbackFcn(app, @GainEditFieldValueChanged, true);
            app.dBmEditField.Tag = 'dBm';
            app.dBmEditField.Placeholder = '0';
            app.dBmEditField.Position = [170 374 45 22];

            % Create GaindBEditFieldLabel
            app.GaindBEditFieldLabel = uilabel(app.UIFigure);
            app.GaindBEditFieldLabel.HorizontalAlignment = 'right';
            app.GaindBEditFieldLabel.Position = [42 336 56 22];
            app.GaindBEditFieldLabel.Text = 'Gain (dB)';

            % Create GaindBEditField
            app.GaindBEditField = uieditfield(app.UIFigure, 'numeric');
            app.GaindBEditField.ValueChangedFcn = createCallbackFcn(app, @GainEditFieldValueChanged, true);
            app.GaindBEditField.Placeholder = '0';
            app.GaindBEditField.Position = [113 336 102 22];

            % Create AttenuationdBEditFieldLabel
            app.AttenuationdBEditFieldLabel = uilabel(app.UIFigure);
            app.AttenuationdBEditFieldLabel.HorizontalAlignment = 'right';
            app.AttenuationdBEditFieldLabel.Position = [6 290 92 22];
            app.AttenuationdBEditFieldLabel.Text = 'Attenuation (dB)';

            % Create AttenuationdBEditField
            app.AttenuationdBEditField = uieditfield(app.UIFigure, 'numeric');
            app.AttenuationdBEditField.ValueChangedFcn = createCallbackFcn(app, @GainEditFieldValueChanged, true);
            app.AttenuationdBEditField.Placeholder = '0';
            app.AttenuationdBEditField.Position = [113 290 102 22];

            % Create MaximumVppEditFieldLabel
            app.MaximumVppEditFieldLabel = uilabel(app.UIFigure);
            app.MaximumVppEditFieldLabel.HorizontalAlignment = 'right';
            app.MaximumVppEditFieldLabel.Position = [344 380 82 22];
            app.MaximumVppEditFieldLabel.Text = 'Maximum Vpp';

            % Create MaximumVppEditField
            app.MaximumVppEditField = uieditfield(app.UIFigure, 'numeric');
            app.MaximumVppEditField.ValueChangedFcn = createCallbackFcn(app, @GainEditFieldValueChanged, true);
            app.MaximumVppEditField.Placeholder = '0';
            app.MaximumVppEditField.Position = [441 380 100 22];
            app.MaximumVppEditField.Value = 2;

            % Create ImpedanceohmsEditField
            app.ImpedanceohmsEditField = uieditfield(app.UIFigure, 'numeric');
            app.ImpedanceohmsEditField.Placeholder = '0';
            app.ImpedanceohmsEditField.Position = [442 339 100 22];
            app.ImpedanceohmsEditField.Value = 50;

            % Create VppEditFieldLabel
            app.VppEditFieldLabel = uilabel(app.UIFigure);
            app.VppEditFieldLabel.HorizontalAlignment = 'right';
            app.VppEditFieldLabel.Position = [122 395 28 22];
            app.VppEditFieldLabel.Text = 'Vpp';

            % Create OutputsLabel
            app.OutputsLabel = uilabel(app.UIFigure);
            app.OutputsLabel.HorizontalAlignment = 'center';
            app.OutputsLabel.FontSize = 24;
            app.OutputsLabel.Position = [114 230 446 32];
            app.OutputsLabel.Text = 'Outputs';

            % Create dBmEditFieldLabel
            app.dBmEditFieldLabel = uilabel(app.UIFigure);
            app.dBmEditFieldLabel.HorizontalAlignment = 'right';
            app.dBmEditFieldLabel.Position = [177 395 31 22];
            app.dBmEditFieldLabel.Text = 'dBm';

            % Create InputsLabel
            app.InputsLabel = uilabel(app.UIFigure);
            app.InputsLabel.HorizontalAlignment = 'center';
            app.InputsLabel.FontSize = 24;
            app.InputsLabel.Position = [114 416 446 32];
            app.InputsLabel.Text = 'Inputs';

            % Create InputLabel
            app.InputLabel = uilabel(app.UIFigure);
            app.InputLabel.Position = [69 374 32 22];
            app.InputLabel.Text = 'Input';

            % Create ImpedanceohmsEditFieldLabel
            app.ImpedanceohmsEditFieldLabel = uilabel(app.UIFigure);
            app.ImpedanceohmsEditFieldLabel.HorizontalAlignment = 'right';
            app.ImpedanceohmsEditFieldLabel.Position = [322 339 105 22];
            app.ImpedanceohmsEditFieldLabel.Text = 'Impedance (ohms)';

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = GainChecker

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