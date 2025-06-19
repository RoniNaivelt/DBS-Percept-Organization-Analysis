function PatientInfo = EditPatientInfo(data, DBS_data)   
    % initialize
    Patiant_Initials = [data.PatientInformation.Initial.PatientFirstName(1) data.PatientInformation.Initial.PatientLastName(1)];
    Patient_name = ['Patient_' Patiant_Initials];
    PatientInfo = DBS_data.(Patient_name).PatientInfo;
end