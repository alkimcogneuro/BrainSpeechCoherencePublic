
% wrote these commands to fix a typo in the EEG_structs, where the audio_file field was pointing to the wrong .wav file.  

load('~/DATA/CodeSwitch_Data/CSepochs/CS_struct06.mat')
EEG_struct.audio_file = "VersionBCriticalAudioParticipant06.wav"
save('~/DATA/CodeSwitch_Data/CSepochs/CS_struct06.mat', 'EEG_struct'); 

load('~/DATA/CodeSwitch_Data/SLepochs/SL_struct06.mat')
EEG_struct.audio_file = "VersionBCriticalAudioParticipant06.wav"
save('~/DATA/CodeSwitch_Data/SLepochs/SL_struct06.mat', 'EEG_struct'); 

load('~/DATA/CodeSwitch_Data/CSepochs/CS_struct10.mat')
EEG_struct.audio_file = "VersionBCriticalAudioParticipant10.wav"
save('~/DATA/CodeSwitch_Data/CSepochs/CS_struct10.mat', 'EEG_struct'); 

load('~/DATA/CodeSwitch_Data/SLepochs/SL_struct10.mat')
EEG_struct.audio_file = "VersionBCriticalAudioParticipant10.wav"
save('~/DATA/CodeSwitch_Data/SLepochs/SL_struct10.mat', 'EEG_struct'); 

load('~/DATA/CodeSwitch_Data/CSepochs/CS_struct20.mat')
EEG_struct.audio_file = "VersionBCriticalAudioParticipant20.wav"
save('~/DATA/CodeSwitch_Data/CSepochs/CS_struct20.mat', 'EEG_struct'); 

load('~/DATA/CodeSwitch_Data/SLepochs/SL_struct20.mat')
EEG_struct.audio_file = "VersionBCriticalAudioParticipant20.wav"
save('~/DATA/CodeSwitch_Data/SLepochs/SL_struct20.mat', 'EEG_struct');  
