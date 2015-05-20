classdef EEGStudio < handle
    
    properties
        sClientDLLPath
        iNch
        fSamplingRate
        viChannelCodes
        caChannelNames
        MainTimer
        mDisplayBuffer
        mProcessingBuffer
        iTimerPeriod
        sTimerCBFunction
        iTimerRepeats
        % State flags
        bDLLLoaded
        bConnected
        bParamsSet 
        bRunning
        iCurrentSample
        hFigureHandle
        hAxisHandle        
        hPlotHandle
        hAxisHandleNFB        
        hPlotHandleNFB
        fYScale
        caChannelTable
        bDisplaySet
        iStartSample
        iStopSample
        UDPObject
        OperationMode % file/device
        DataFileHandle
        DataFileName
        FileData
        FilePos
        FileSize
    end % properties
     
    methods
%        function eegs = EEGStudio()
%        function SetDLLPath(eegs, sPath)
%        function LoadDLL(eegs)
%        function Connect(eegs)
%        function Disconnect(eegs)
%        function GetAcqParams(eegs)
%        function InitTimer(eegs)
%        function Start(eegs)
%        function Stop(eegs)
%        function LoadDLLAndConnect(eegs, sDLLPath)
%        function MasterStart(eegs, sDLLPath)
%        function MasterStop(eegs)
%        function InitDisplay(eegs, hUserFigureHandle)
%        function StopDisplay(eegs)
%        function PlotEEGData(eegs, nstart, n2plot)
%        function NewData = GetNewData(eegs) 
%        function SaveSessionData(eegs,sPath)
%        function DefaultNewData(eegs, event)
%        funciton InitUDP(eegs,rHost, rPort);
        
        % CONSTRUCTOR
        function eegs = EEGStudio()
            eegs.sClientDLLPath  = '';
            eegs.iNch = -1;
            eegs.fSamplingRate = -1.0;
            eegs.viChannelCodes = [];
            eegs.caChannelNames = {};
            eegs.MainTimer = [];
            eegs.mDisplayBuffer = [];
            eegs.mProcessingBuffer =[]; 
            eegs.iTimerPeriod = 0.1; % default to 0.1 second
            eegs.sTimerCBFunction = 'DefaultNewData';
            eegs.iTimerRepeats = Inf;
            eegs.bConnected = false;
            eegs.bDLLLoaded = false;
            eegs.bParamsSet = false;
            eegs.bRunning   = false;
            eegs.iCurrentSample = 0;
            eegs.hFigureHandle  = []; 
            eegs.hAxisHandle = [];
            eegs.hPlotHandle = [];
            eegs.hPlotHandleNFB = [];
            eegs.hAxisHandleNFB = [];
            eegs.fYScale = 3000;
            eegs.bDisplaySet = false;
            eegs.caChannelTable = InitChannelTable();
            eegs.iStartSample = 0;
            eegs.iStopSample = 0;
            eegs.UDPObject = [];
            eegs.OperationMode =   1; % device
            eegs.DataFileHandle = -1;
            eegs.DataFileName = '';
        end
        
        function InitUDP(eegs,rHost, rPort,lPort)
            a = instrfindall;
            if(~isempty(a))
                fclose(a);
            end;
           if(nargin==1) 
               rPort = 9090;
               lPort = 9091;
               rHost = '192.168.2.2';
           end;
           eegs.UDPObject = udp(rHost,'RemotePort',rPort,'LocalPort', lPort);
           fopen(eegs.UDPObject);
        end
 
        function StopUDP(eegs)
           echoudp('off');
           fclose(eeg.UDPObject);
        end

        function SetDLLPath(eegs, sPath)
            % check if the path is valid
            bOK = exist(sPath,'file');
            if(bOK==0)
                MException('EEGStudio:SetDLLPath', ...
                           'InvalidPath %s.', sPath);
                return;
            end
            % if path is valid - assign
            eegs.sClientDLLPath = sPath;
        end;
        
        function LoadDLL(eegs)
            if(isempty(eegs.sClientDLLPath))
                MException('EEGStudio:LoadDLL', ...
                           'No DLL path set. Use ::SetDLLPath first.');
            end;
           
            try
                NET.addAssembly(eegs.sClientDLLPath);
            catch ME
                MException(ME);
            end;
            
            eegs.bDLLLoaded = true;
        end;
        
        function Connect(eegs)
            
            if(eegs.OperationMode==1) % device
                if(~eegs.bDLLLoaded)
                    MException('EEGStudio:Connect', ...
                               'No DLL loaded. Use ::LoadDLL first.');
                end;
                EEGStudio_Client.Methods.Start;
                eegs.bConnected = true;
            else
                % look for file
                if(exist(eegs.DataFileName,'file'))
                    eegs.DataFileHandle = fopen(eegs.DataFileName,'rb');
                    if(eegs.DataFileHandle > 0)
                        eegs.bConnected = true;
                    else
                        eegs.bConnected = false;
                    end;
                else
                    MException('EEGStudio:Connect', ...
                               'FileMode: no datafile found');
                end
            end
            
        end;
        
        function Disconnect(eegs)
            if(eegs.bConnected)
                EEGStudio_Client.Methods.Stop;
                eegs.bConnected = false;
            end;
        end;
        
        function GetAcqParams(eegs)
            if(eegs.bConnected)
                if(eegs.OperationMode==1) %device
                    eegs.iNch = EEGStudio_Client.Methods.GetChannelsCount();
                    if(eegs.iNch==0)
                        eegs.bParamsSet= false;
                        return;
                    end;
                    eegs.fSamplingRate = EEGStudio_Client.Methods.GetSamplingFrequency();
                    eegs.viChannelCodes = zeros(1,eegs.iNch);
                    eegs.caChannelNames = cell(1,eegs.iNch);
                    for i=1:eegs.iNch
                        eegs.viChannelCodes(i) = EEGStudio_Client.Methods.GetChannelCode(i-1);
                        eegs.caChannelNames{i} = eegs.caChannelTable{eegs.viChannelCodes(i)+1}; 
                        % TODO Create UIDS table and set names according to it
                    end;
                    eegs.bParamsSet= true;
                else % if file
                    frewind(eegs.DataFileHandle);
                    eegs.iNch = fread(eegs.DataFileHandle,1,'int');
                    eegs.fSamplingRate  = fread(eegs.DataFileHandle,1,'float');
                    eegs.viChannelCodes = zeros(1,eegs.iNch);
                    eegs.caChannelNames = cell(1,eegs.iNch);
                    for i=1:eegs.iNch
                        eegs.caChannelNames{i} = num2str(i); 
                    end;
                    eegs.bParamsSet= true;
                end;
            end;
        end;
        
        function InitTimer(eegs, sMyCBFunction)
            % if a timer object is present start delete the old one and start anew
            if(isa(eegs.MainTimer,'timer'))
                stop(eegs.MainTimer);
                delete(eegs.MainTimer);
            end;
            if(nargin==1)
                eegs.MainTimer = timer('TimerFcn',@(src, evt)eegs.DefaultNewData);
            else
                eegs.MainTimer = timer('TimerFcn',sMyCBFunction);
            end;
            set(eegs.MainTimer,'ExecutionMode','fixedSpacing','Period',eegs.iTimerPeriod, 'TasksToExecute',eegs.iTimerRepeats);
        end;
                
        function Start(eegs)
            if(~isa(eegs.MainTimer,'timer'))
                MException('EEGStudio:Start', ...
                                'MainTimer is not initialized. Use ::InitTimer first.');
            end;
 
            if(~eegs.bConnected)
                MException('EEGStudio:Start', ...
                                'Not connected. Use ::Connect first.');
            end;
            
            if(~eegs.bParamsSet)
                MException('EEGStudio:Start', ...
                                'ACQ params not set. Use ::GetAcaParams first.');
            end;
            
            if(eegs.OperationMode==0)
                % load data from file into memory
                frewind(eegs.DataFileHandle);
                % skip two header params
                fread(eegs.DataFileHandle,1,'int');
                fread(eegs.DataFileHandle,1,'float');
                k=1;
                try
                    while(~feof(eegs.DataFileHandle))
                        eegs.FileData(:,k)= fread(eegs.DataFileHandle,eegs.iNch,'real*4');
                        k=k+1;
                    end;
                catch
                end;
                fclose(eegs.DataFileHandle);
                eegs.iStartSample = 1;
                eegs.FilePos = 1;
                eegs.FileSize = size(eegs.FileData,2);
            else
                eegs.iStartSample = EEGStudio_Client.Methods.GetTicksCount();
            end;
            
            start(eegs.MainTimer);
            eegs.bRunning = true;
        end;    

        function Stop(eegs)
            if(~eegs.bRunning)
                return;
            end;
            if(~isa(eegs.MainTimer,'timer'))
                MException('EEGStudio:Start', ...
                                'MainTimer is not initialized. Use ::InitTimer first.');
            end;
            
            stop(eegs.MainTimer);
            eegs.iStopSample = EEGStudio_Client.Methods.GetTicksCount();
            nss.bRunning = false;
        end;
        
        function SetDataFileName(eegs, sDataFilePath)
            eegs.DataFileName = sDataFilePath;
        end;
        
        function LoadDLLAndConnect(eegs, sDLLPath)
            if(eegs.OperationMode==1)
                eegs.SetDLLPath(sDLLPath);
                eegs.LoadDLL();
            else
                eegs.DataFileName = sDLLPath;
            end;
            eegs.Connect();
        end;
            
        function MasterStart(eegs, sDLLPath)

            eegs.LoadDLLAndConnect(sDLLPath)
            % Now we will obtain ACQ parameters and set the corresponding
            % fields of this object used later in the timer callback
            % functions
            eegs.GetAcqParams();
            %You may want to change main timer settings 
            % iTimerPeriod, iTimerCBFunction, iTimerRepeats here before calling
            %InitTimer function
            eegs.InitTimer();
            % now we are ready to start
            eegs.Start();
            eegs.InitDisplay();
            disp('Locate and press MatLab button under the Add-on tab of EEG studio to start data transfer.');
        end;
        
        function MasterStartNFB(eegs, sDLLPath)

            eegs.LoadDLLAndConnect(sDLLPath)
            % Now we will obtain ACQ parameters and set the corresponding
            % fields of this object used later in the timer callback
            % functions
            eegs.GetAcqParams();
            %You may want to change main timer settings 
            % iTimerPeriod, iTimerCBFunction, iTimerRepeats here before calling
            %InitTimer function
            eegs.iTimerPeriod = 0.1;
            eegs.sTimerCBFunction = 'eegs.NewDataNFB';
            eegs.InitTimer(eegs.sTimerCBFunction);
            % now we are ready to start
            eegs.Start();
            eegs.InitDisplayNFB();
            disp('Locate and press MatLab button under the Add-on tab of EEG studio to start data transfer.');
        end;
        
        function MasterStartUDP(eegs, sDLLPath)

            eegs.LoadDLLAndConnect(sDLLPath)
            % Now we will obtain ACQ parameters and set the corresponding
            % fields of this object used later in the timer callback
            % functions
            eegs.GetAcqParams();
            %You may want to change main timer settings 
            % iTimerPeriod, iTimerCBFunction, iTimerRepeats here before calling
            %InitTimer function
            eegs.iTimerPeriod = 0.1;
            eegs.sTimerCBFunction = 'eegs.NewDataUDP';
            eegs.InitTimer(eegs.sTimerCBFunction);
            % now we are ready to start
            eegs.InitUDP();
            eegs.Start();
            eegs.InitDisplayNFB(); % use the same display
            disp('Locate and press MatLab button under the Add-on tab of EEG studio to start data transfer.');
        end;
        
        function MasterStop(eegs)
            eegs.StopDisplay;
            eegs.Stop;
            eegs.Disconnect();
            eegs.iCurrentSample = 0;
            if(~isempty(eegs.UDPObject))
                eegs.StopUDP();
            end
%            eegs.del;
        end;
        % DESTRUCTOR.
        function del(eegs)
            eegs.Stop();
            eegs.Disconnect();
            delete(eegs);
        end  
        
        function ConfigureEEGAxis(eegs, hAxis,n2plot)
            ax = hAxis; % get(eegs.hFigureHandle,'Children');
            set(ax,'YLimMode', 'manual');
            set(ax,'YLim',[0 eegs.iNch+1]*eegs.fYScale);
            set(ax,'XLimMode', 'manual');
            set(ax,'XLim',[1 n2plot]);
            set(ax,'Yticklabel',[]);
            set(ax,'Xticklabel',[]);
            set(ax,'Ytick', (( 1:eegs.iNch ))*eegs.fYScale );
            set(ax,'Yticklabel',eegs.caChannelNames);
            set(ax,'YDir','reverse'); % first channel on the top
            set(ax,'DrawMode','fast');
            set(ax,'Xtick',1:50:n2plot);
            XTick = get(ax,'XTick');
            XTickLabel = cell(1,length(XTick));
            for i=1:length(XTick)
                XTickLabel{i} = sprintf('%2.2f',(single(XTick(i)))/single(eegs.fSamplingRate));
            end;
            set(ax,'XtickLabel',XTickLabel);
            set(ax,'FontSize',8);
            xlabel('time, ms');
            title('EEG data stream');
        end
        
        function InitDisplay(eegs, hUserFigureHandle)
            
            if(eegs.bDisplaySet)
                return;
            end;
            if(~eegs.bParamsSet)
                MException('EEGStudio:InitDisplay', ...
                                'ACQ params not set. Use ::GetAcqParams first.');
            end;
            
            if(eegs.iNch==0)
                eegs.bDisplaySet = false;
                return;
            end;
            
            if(nargin==2)
                eegs.hFigureHandle = hUserFigureHandle;
            else
                eegs.hFigureHandle = figure;
            end;
            
            n2plot = 10*fix(eegs.iTimerPeriod*eegs.fSamplingRate);
            eegs.mDisplayBuffer = zeros(eegs.iNch, n2plot);
            eegs.hAxisHandle = subplot(1,1,1,'Parent',eegs.hFigureHandle);
            eegs.hPlotHandle = plot(eegs.mDisplayBuffer','Parent',eegs.hAxisHandle);
            eegs.ConfigureEEGAxis(eegs.hAxisHandle,n2plot);
            eegs.bDisplaySet = true;
        end
        
        function InitDisplayNFB(eegs, hUserFigureHandle)
            
            if(eegs.bDisplaySet)
                return;
            end;
            if(~eegs.bParamsSet)
                MException('EEGStudio:InitDisplay', ...
                                'ACQ params not set. Use ::GetAcqParams first.');
            end;
            
            if(eegs.iNch==0)
                eegs.bDisplaySet = false;
                return;
            end;
            
            if(nargin==2)
                eegs.hFigureHandle = hUserFigureHandle;
            else
                eegs.hFigureHandle = figure;
            end;
            
            n2plot = 10*fix(eegs.iTimerPeriod*eegs.fSamplingRate);
            eegs.mDisplayBuffer = zeros(eegs.iNch, n2plot);
            eegs.hAxisHandle = subplot(1,2,1,'Parent',eegs.hFigureHandle);
            eegs.hPlotHandle = plot(eegs.mDisplayBuffer','Parent',eegs.hAxisHandle);
            eegs.ConfigureEEGAxis(eegs.hAxisHandle,n2plot);
            
            eegs.hAxisHandleNFB = subplot(1,2,2,'Parent',eegs.hFigureHandle);
            eegs.hPlotHandleNFB = bar(1,1,'Parent',eegs.hAxisHandleNFB);
            set(eegs.hAxisHandleNFB,'XLim',[0 2]);
            set(eegs.hAxisHandleNFB,'YLim',[0 10]);
            set(eegs.hAxisHandleNFB,'FontSize',8);
            title('Feedback sgnal = alpha/theta @ P4');
            eegs.bDisplaySet = true;
        end

        
        function StopDisplay(eegs)
            try
                if(~isempty(eegs.hAxisHandle))
                    if(~isempty(eegs.hFigureHandle))
                        close(eegs.hFigureHandle);
                        eegs.hFigureHandle = [];
                    end;
                    eegs.hAxisHandle   = [];
                    eegs.hPlotHandle   = [];
                end;
            catch
                %no catch actions
            end;
        end;
        
        function PlotEEGData(eegs, nstart, n2plot)
            if(length(eegs.hPlotHandle)~=eegs.iNch)
                return;
            end;
            
            if(eegs.OperationMode==1)
                for i=1:eegs.iNch
                    set(eegs.hPlotHandle(i),'YData',...
                    fix(eegs.fYScale*double(i)) + single(EEGStudio_Client.Methods.GetRange(i-1, nstart, n2plot)));                    
                    set(eegs.hPlotHandle(i),'XData', 1:n2plot);                    
                end;
            else %use data loaded from file in Start() function
                if(nstart<=0)
                    Data =  eegs.FileData(:,eegs.FileSize+nstart+1:eegs.FileSize);
                    Data = [Data eegs.FileData(:,1:n2plot+nstart)];
                else    
                    Data = eegs.FileData(:,eegs.FilePos-n2plot+1:eegs.FilePos);
                end

                for i=1:eegs.iNch
                    set(eegs.hPlotHandle(i),'YData',...
                    fix(eegs.fYScale*double(i)) + 5*Data(i,:));                    
                    set(eegs.hPlotHandle(i),'XData', 1:n2plot);                    
                end;
            end;
            
        end;
        
        function SaveSessionData(eegs,sPath)

            if nargin~=2
                disp('SaveSessionData: Incorrect number of arguments');
                return;
            end;
            
            n2save = eegs.iStopSample - eegs.iStartSample;
            
            Data = zeros(eegs.iNch, n2save);
            
            for i=1:eegs.iNch
                Data(i,1:n2save) = single(EEGStudio_Client.Methods.GetRange(i-1, eegs.iStartSample, n2save));                    
            end;
            
            save(sPath,'Data','eegs');
            
        end;    
        
        function Data = GetSessionData(eegs)

            if nargin~=2
                disp('SaveSessionData: Incorrect number of arguments');
                return;
            end;
            
            n2get = eegs.iStopSample - eegs.iStartSample;
            
            Data = zeros(eegs.iNch, n2get);
            
            for i=1:eegs.iNch
                Data(i,1:n2save) = single(EEGStudio_Client.Methods.GetRange(i-1, eegs.iStartSample, n2get));                    
            end;
            
        end;    
            
        function DefaultNewData(eegs, event)
            %In order to ensure good synchronization we recomment that
            %Server is started past the moment when the client has been
            %started. Therefore, it is possible that parameters are not
            %set.
            
            if(~eegs.bParamsSet)
                eegs.GetAcqParams();
                eegs.InitDisplay();
            end;
            if(eegs.OperationMode==1)
                nav = EEGStudio_Client.Methods.GetTicksCount();
                %%this is a workaround the existing problem
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if(nav<eegs.iStartSample)
                    eegs.iStartSample = 0;
                end;
                n2plot = 10*fix(eegs.iTimerPeriod*eegs.fSamplingRate);
                if(length(eegs.hPlotHandle)==eegs.iNch && nav>n2plot)
                    eegs.PlotEEGData(nav-n2plot, n2plot);
                end;
                % update current sample value only once and in the end
                eegs.iCurrentSample = nav;
            else
                if(length(eegs.hPlotHandle)==eegs.iNch)
                    if(eegs.FilePos>n2plot)
                        eegs.PlotEEGData(eegs.FilePos-n2plot,n2plot);
                    end;
                    eegs.FilePos = mod(eegs.FilePos+fix(eegs.iTimerPeriod*eegs.fSamplingRate), eegs.FileSize);
                end;
                eegs.iCurrentSample = eegs.iCurrentSample + fix(eegs.iTimerPeriod*eegs.fSamplingRate);
            end;
                
        end;
        
        function NewDataNFB(eegs, event)
            %In order to ensure good synchronization we recomment that
            %Server is started past the moment when the client has been
            %started. Therefore, it is possible that parameters are not
            %set.
            if(~eegs.bParamsSet)
                eegs.GetAcqParams();
                eegs.InitDisplayNFB();
            end;
            
            if(eegs.OperationMode==1)
                nav = EEGStudio_Client.Methods.GetTicksCount();
                if(nav<eegs.iStartSample)
                    eegs.iStartSample = 0;
                end;
            end;
            
            %%this is a workaround the existing problem
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            n2plot = 10*fix(eegs.iTimerPeriod*eegs.fSamplingRate);
            n2analyze = 2*fix(eegs.iTimerPeriod*eegs.fSamplingRate);            
            iChannelIndexNFB =14; % P4 
            if(eegs.OperationMode==1)
                if(length(eegs.hPlotHandle)==eegs.iNch && nav>n2plot)
                    eegs.PlotEEGData(nav-n2plot, n2plot);
                end;
                Data = single(EEGStudio_Client.Methods.GetRange(iChannelIndexNFB,nav-n2analyze,n2analyze));
                % updatecurrent sample value only once and in the end
                eegs.iCurrentSample = nav;
            else
                if(length(eegs.hPlotHandle)==eegs.iNch)
                    if(eegs.FilePos>n2plot)
                        eegs.PlotEEGData(eegs.FilePos-n2plot,n2plot);
                    end;
                end;
                if(eegs.FilePos-n2analyze+1<=0)
                   n = abs(eegs.FilePos-n2analyze);
                   Data =  eegs.FileData(:,eegs.FileSize-n+1:eegs.FileSize);
                   Data = [Data eegs.FileData(:,1:n2analyze-n)];
                else    
                    Data = eegs.FileData(:,eegs.FilePos-n2analyze+1,eegs.FilePos);
                end
                eegs.FilePos = mod(eegs.FilePos+fix(eegs.iTimerPeriod*eegs.fSamplingRate), eegs.FileSize);
                eegs.iCurrentSample = eegs.iCurrentSample + fix(eegs.iTimerPeriod*eegs.fSamplingRate);
            end;
            
            DataFFT = fft(Data);
            %The following 3 lines don't have to be executed every callback funciton
            % and can be precomputed
            freq = single(eegs.fSamplingRate)*single(1:length(Data)-1)/(length(Data)-1);
            alpha_range = find(freq>=8 & freq<=12 );
            theta_range = find(freq>=4 & freq<8 );
            Q = mean(abs(DataFFT(alpha_range)))./mean(abs(DataFFT(theta_range)));
            set(eegs.hPlotHandleNFB,'YData',Q);
        end;
        
        function NewDataUDP(eegs, event)
            %In order to ensure good synchronization we recomment that
            %Server is started past the moment when the client has been
            %started. Therefore, it is possible that parameters are not
            %set.
            if(~eegs.bParamsSet)
                eegs.GetAcqParams();
                eegs.InitDisplayNFB();
            end;
            
            if(eegs.OperationMode==1)
                nav = EEGStudio_Client.Methods.GetTicksCount();
                if(nav<eegs.iStartSample)
                    eegs.iStartSample = 0;
                end;
            end;
            
            %%this is a workaround the existing problem
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            n2plot = 10*fix(eegs.iTimerPeriod*eegs.fSamplingRate);
            n2analyze = 2*fix(eegs.iTimerPeriod*eegs.fSamplingRate);            
            iChannelIndexNFB =14; % P4 
            if(eegs.OperationMode==1)
                if(length(eegs.hPlotHandle)==eegs.iNch && nav>n2plot)
                    eegs.PlotEEGData(nav-n2plot, n2plot);
                end;
                Data = single(EEGStudio_Client.Methods.GetRange(iChannelIndexNFB,nav-n2analyze,n2analyze));
                % updatecurrent sample value only once and in the end
                eegs.iCurrentSample = nav;
            else
                if(length(eegs.hPlotHandle)==eegs.iNch)
                    if(eegs.FilePos>n2plot)
                        eegs.PlotEEGData(eegs.FilePos-n2plot,n2plot);
                    end;
                end;
                if(eegs.FilePos-n2analyze+1<=0)
                   n = abs(eegs.FilePos-n2analyze);
                   Data =  eegs.FileData(:,eegs.FileSize-n+1:eegs.FileSize);
                   Data = [Data eegs.FileData(:,1:n2analyze-n)];
                else    
                    Data = eegs.FileData(:,eegs.FilePos-n2analyze+1:eegs.FilePos);
                end
                eegs.FilePos = mod(eegs.FilePos+fix(eegs.iTimerPeriod*eegs.fSamplingRate), eegs.FileSize)+1;
                eegs.iCurrentSample = eegs.iCurrentSample + fix(eegs.iTimerPeriod*eegs.fSamplingRate);
            end;
            
            DataFFT = fft(Data);
            %The following 3 lines don't have to be executed every callback funciton
            % and can be precomputed
            freq = single(eegs.fSamplingRate)*single(1:length(Data)-1)/(length(Data)-1);
            alpha_range = find(freq>=8 & freq<=12 );
            theta_range = find(freq>=4 & freq<8 );
            Q = mean(abs(DataFFT(alpha_range)))./mean(abs(DataFFT(theta_range)));
            set(eegs.hPlotHandleNFB,'YData',Q);
            if(~isempty(eegs.UDPObject))
                fwrite(eegs.UDPObject,repmat(Q,1,100));
            end;
        end;
        
       
        % get the data aquired between two consequent callback calls
        function NewData = GetNewData(eegs) 
            nav = EEGStudio_Client.Methods.GetTicksCount();
            n2get = nav-eegs.iCurrentSample;
            NewData = zeros(eegs.iNch, n2get);
            for i=1:eegs.iNch
                NewData(i,1:n2get) = single(EEGStudio_Client.Methods.GetRange(i-1, eegs.iCurrentSample, n2get));                    
            end;
        end;
        
        
    end  % methods   
end
 

