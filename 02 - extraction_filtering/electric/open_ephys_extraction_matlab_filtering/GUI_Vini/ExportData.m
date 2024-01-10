function ExportData(DATA, channels, directory, FileName)
%function ExportData(DATA, channels, directory)
%
%Exports and saves data recorded by Open Ephys in Binary format
%  DATA: Data struucture read by load_open_ephys_bynary 
%  channels: Array of channels to save.
%  directory: The directory to save the file
%  FileName: Name of the file to save
%
%Example:
% DATA = DATA
% channels = [1:10, 15] % To save the channels 1 to 10 and 15
% Directory = C:\User\documents\experiments\electric
% FileName = ElectricData
%


            app.channels = str2num(app.ChannelsEditField.Value); 
            FileName = fullfile(app.ExpPath, strcat(app.FileNameEditField_2.Value, '.mat'));

            if app.RawDataButton.Value == 1
                expDataRaw.Header.Channels = app.channels;
                expDataRaw.Header.Fs = app.Fs;
                expDataRaw.Header.Recording = app.RecordingEditField.Value;
                expDataRaw.Header.Date = app.DateEditField.Value;
                expDataRaw.Header.OBS = app.OBSTextArea.Value;
                expDataRaw.Header.Animal = app.AnimalEditField.Value;
                expDataRaw.Header.Mass = app.MassEditField.Value;

                expDataRaw.Data = app.OData(app.channels, :);
                expDataRaw.Timestamps = app.D.Timestamps;
                expDataRaw.TTL = app.E.Timestamps(1:2:end);
                
                save(FileName, '-struct', 'expDataRaw');

            elseif app.FilteredDataButton.Value == 1 
                expDataFil.Header.Channels = app.channels;
                expDataFil.Header.Fs = app.Fs;
                expDataFil.Header.Recording = app.RecordingEditField.Value;
                expDataFil.Header.Date = app.DateEditField.Value;
                expDataFil.Header.OBS = app.OBSTextArea.Value;
                expDataFil.Header.Animal = app.AnimalEditField.Value;
                expDataFil.Header.Mass = app.MassEditField.Value;

                expDataFil.Data = app.OData(app.channels, :);
                expDataFil.Timestamps = app.D.Timestamps;
                expDataFil.TTL = app.E.Timestamps(1:2:end);

                har_Q = app.QualityEditField.Value;
                har_freq = app.NaturalFrequencyEditField.Value;
                har_order = app.OrderEditField.Value;
                hp_freq = app.CutoffFrequencyEditField_2.Value;
                lp_freq = app.CutoffFrequencyEditField.Value;
                sg_order = app.OrderEditField_2.Value;
                sg_framelen = app.FrameLengthoddEditField.Value;
    
                checkboxValue1_4 = [app.ActiveCheckBox.Value
                    app.ActiveCheckBox_2.Value
                    app.ActiveCheckBox_3.Value
                    app.ActiveCheckBox_4.Value];               

                for i = 1:length(app.channels)
                    expDataFil.Data(i, :) = filterdata(expDataFil.Data(i, :), checkboxValue1_4, har_Q, har_freq, har_order,...
                        hp_freq, lp_freq, sg_order, sg_framelen, app.Fs);
                end

                save(FileName, '-struct', 'expDataFil'); 
            end                           
        end