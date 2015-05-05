classdef SubjectRecord < handle
    properties
        subject_name
        acq_date
        time_start
        time_stop
        folder_name
        nss_file_name
        summary_file_name
        this_file_name
    end 
   
    methods
        function self = SubjectRecord(self) %#ok<INUSD>
            self.subject_name       = 'Null';
            self.acq_date           = date;
            self.time_start         = datestr(now,'HH-MM-SS');
            self.time_stop          = 0;
            self.folder_name        = '';
            self.nss_file_name      = '';
            self.summary_file_name  = '';
            self.this_file_name     = '';
        end
    end
end


