#upload some libraries
import os
import sys
import subprocess
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib import style
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from PIL import Image
import datetime
import math
from scipy import stats
import statsmodels.formula.api as smf
import statsmodels.api as sm
import statistics
import glob
from os.path import splitext, basename
import sqlite3

st.set_page_config(
    layout = 'wide'
)
# set paths
twasdir = 'T-WAS'
pwasdir = 'P-WAS'
tewasdir = 'TE-WAS'

#read omics data files
#chromosome list
#chr_list = ['Chr' + str(i) for i in range(1, 23)]
#num_list = [i for i in range(1, 23)]
chr_list = [str(i) for i in range(1, 23)]

#twas_case
twas_case_folder = f'{twasdir}/Top_case/'
twas_case_file = glob.glob(twas_case_folder + "PD.*.top")
twas_case = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in twas_case_file}

#twas_control
twas_control_folder = f'{twasdir}/Top_control/'
twas_control_file = glob.glob(twas_control_folder + "PD_control.*.top")
twas_control = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in twas_control_file}


#tewas_case
tewas_case_folder = f'{tewasdir}/Top_case/'
tewas_case_file = glob.glob(tewas_case_folder + "PD.*.top")
tewas_case = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in tewas_case_file}

#tewas_control
tewas_control_folder = f'{tewasdir}/Top_control/'
tewas_control_file = glob.glob(tewas_control_folder + "PD_control.*.top")
tewas_control = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in tewas_control_file}


#pwas_csf_cardio_case
pwas_csf_cardio_case_folder = f'{pwasdir}/Top_CSF_Cardio_case/'
pwas_csf_cardio_case_file = glob.glob(pwas_csf_cardio_case_folder + "PD.*.top")
pwas_csf_cardio_case = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_csf_cardio_case_file}

#pwas_csf_cardio_control
pwas_csf_cardio_control_folder = f'{pwasdir}/Top_CSF_Cardio_control/'
pwas_csf_cardio_control_file = glob.glob(pwas_csf_cardio_control_folder + "PD.*.top")
pwas_csf_cardio_control = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_csf_cardio_control_file}

#pwas_csf_inf_case
pwas_csf_inf_case_folder = f'{pwasdir}/Top_CSF_INF_case/'
pwas_csf_inf_case_file = glob.glob(pwas_csf_inf_case_folder + "PD.*.top")
pwas_csf_inf_case = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_csf_inf_case_file}

#pwas_csf_inf_control
pwas_csf_inf_control_folder = f'{pwasdir}/Top_CSF_INF_control/'
pwas_csf_inf_control_file = glob.glob(pwas_csf_inf_control_folder + "PD.*.top")
pwas_csf_inf_control = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_csf_inf_control_file}

#pwas_csf_neuro_case
pwas_csf_neuro_case_folder = f'{pwasdir}/Top_CSF_Neuro_case/'
pwas_csf_neuro_case_file = glob.glob(pwas_csf_neuro_case_folder + "PD.*.top")
pwas_csf_neuro_case = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_csf_neuro_case_file}

#pwas_csf_neuro_control
pwas_csf_neuro_control_folder = f'{pwasdir}/Top_CSF_Neuro_control/'
pwas_csf_neuro_control_file = glob.glob(pwas_csf_neuro_control_folder + "PD.*.top")
pwas_csf_neuro_control = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_csf_neuro_control_file}

#pwas_csf_onc_case
pwas_csf_onc_case_folder = f'{pwasdir}/Top_CSF_ONC_case/'
pwas_csf_onc_case_file = glob.glob(pwas_csf_onc_case_folder + "PD.*.top")
pwas_csf_onc_case = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_csf_onc_case_file}

#pwas_csf_onc_control
pwas_csf_onc_control_folder = f'{pwasdir}/Top_CSF_ONC_control/'
pwas_csf_onc_control_file = glob.glob(pwas_csf_onc_control_folder + "PD.*.top")
pwas_csf_onc_control = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_csf_onc_control_file}


#pwas_plasma_cardio_case
pwas_plasma_cardio_case_folder = f'{pwasdir}/Top_Plasma_Cardio_case/'
pwas_plasma_cardio_case_file = glob.glob(pwas_plasma_cardio_case_folder + "PD.*.top")
pwas_plasma_cardio_case = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_plasma_cardio_case_file}

#pwas_plasma_cardio_control
pwas_plasma_cardio_control_folder = f'{pwasdir}/Top_Plasma_Cardio_control/'
pwas_plasma_cardio_control_file = glob.glob(pwas_plasma_cardio_control_folder + "PD.*.top")
pwas_plasma_cardio_control = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_plasma_cardio_control_file}

#pwas_plasma_inf_case
pwas_plasma_inf_case_folder = f'{pwasdir}/Top_Plasma_INF_case/'
pwas_plasma_inf_case_file = glob.glob(pwas_plasma_inf_case_folder + "PD.*.top")
pwas_plasma_inf_case = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_plasma_inf_case_file}

#pwas_plasma_inf_control
pwas_plasma_inf_control_folder = f'{pwasdir}/Top_Plasma_INF_control/'
pwas_plasma_inf_control_file = glob.glob(pwas_plasma_inf_control_folder + "PD.*.top")
pwas_plasma_inf_control = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_plasma_inf_control_file}

#pwas_plasma_neuro_case
pwas_plasma_neuro_case_folder = f'{pwasdir}/Top_Plasma_Neuro_case/'
pwas_plasma_neuro_case_file = glob.glob(pwas_plasma_neuro_case_folder + "PD.*.top")
pwas_plasma_neuro_case = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_plasma_neuro_case_file}

#pwas_plasma_neuro_control
pwas_plasma_neuro_control_folder = f'{pwasdir}/Top_Plasma_Neuro_control/'
pwas_plasma_neuro_control_file = glob.glob(pwas_plasma_neuro_control_folder + "PD.*.top")
pwas_plasma_neuro_control = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_plasma_neuro_control_file}

#pwas_plasma_onc_case
pwas_plasma_onc_case_folder = f'{pwasdir}/Top_Plasma_ONC_case/'
pwas_plasma_onc_case_file = glob.glob(pwas_plasma_onc_case_folder + "PD.*.top")
pwas_plasma_onc_case = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_plasma_onc_case_file}

#pwas_plasma_onc_control
pwas_plasma_onc_control_folder = f'{pwasdir}/Top_Plasma_ONC_control/'
pwas_plasma_onc_control_file = glob.glob(pwas_plasma_onc_control_folder + "PD.*.top")
pwas_plasma_onc_control = {splitext(basename(file))[0].split('.')[-1] : pd.read_csv(file, sep='\t', index_col=False, usecols =['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P']) for file in pwas_plasma_onc_control_file}

#empty dataframe 
df_1 = pd.DataFrame(columns = ['ID', 'CHR', 'EQTL.ID', 'EQTL.Z', 'TWAS.Z', 'TWAS.P'])

#data for ternary plot
df = pd.read_csv('ternary_plot.csv')

##Background color
def local_css(file_name):
    with open(file_name) as f:
        st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
local_css("style.css")
####################### HEAD ##############################################

head_1, title, head_4 = st.beta_columns([1.2, 4, 1])

gp2 = Image.open('gp2_2.jpg')
head_1.image(gp2, width = 120)

with title:
    st.markdown("""
    <style>
    .big-font {
        font-family:IBM Plex Sans; color:#0f557a; font-size:48px !important;
    }
    </style>
    """, unsafe_allow_html=True)
    st.markdown('<p class="big-font">Multi-omics Results</p>', unsafe_allow_html=True)


with head_4:
    def modification_date(filename):
        t = os.path.getmtime(filename)
        return datetime.datetime.fromtimestamp(t)
    date = modification_date("ternary_plot.csv") 
    st.markdown("**RESULT DATE**")  
    st.markdown(date)

    
    
#DB Management
conn = sqlite3.connect('data.db')
c = conn.cursor()

def create_usertable():
    c.execute('CREATE TABLE IF NOT EXISTS usertable(username TEXT, password TEXT)')
    
def add_userdata(username, password):
    c.execute('INSERT INTO usertable(username, password) VALUES (?,?)', (username, password))
    conn.commit()
    
def login_user(username, password):
    c.execute('SELECT * FROM usertable WHERE username =? AND password = ?', (username, password))
    data = c.fetchall()
    return data
    
def view_all_users():
    c.execute('SELECT * FROM usertable')
    data = c.fetchall()
    return data




def main():
    """"Simple Login App"""
    st.title("Simple Login App")
    
    menu = ["Home", "Login", "SignUp"]
    choice = st.sidebar.selectbox("Menu", menu)
    
    if choice == "Home":
        st.subheader("Home")
        
    elif choice == "Login":
        st.subheader("Login Section")
        
        username = st.sidebar.text_input("User Name")
        password = st.sidebar.text_input("Password", type='password')
        if st.sidebar.checkbox("Login"):
            #if password == '12345':
            create_usertable()
            result = login_user(username, password)
            if result: 
                st.success("Logged in as {}".format(username))       
    
    

                ########################  SIDE BAR #########################################
                st.sidebar.markdown('**Data Selection**', unsafe_allow_html=True)
                #prune_selection = df_qc['step'].unique().tolist()
                #prune_selection.remove('variant_prune')
                selected_metrics = st.sidebar.selectbox(label="T-WAS Chomosome Number", options=chr_list)
                #selected_metrics = df_qc.loc[df_qc['step'] == prune_selection]
                #selected_metrics = selected_metrics.reset_index()
                selected_metrics_1 = st.sidebar.selectbox(label="TE-WAS Chomosome Number", options=chr_list)

                selected_metrics_2 = st.sidebar.selectbox(label="P-WAS Sample", options=['CSF', 'Plasma'])
                selected_metrics_3 = st.sidebar.selectbox(label="P-WAS Panel", options=['Cardiometabolic',  'Inflammation', 'Neurology', 'Oncology'])
                selected_metrics_4 = st.sidebar.selectbox(label="P-WAS Chomosome Number", options=chr_list)

                ########################  right column  #########################################
                left_column, right_column = st.beta_columns([1.2,1.2])
                with right_column: 
                    #twas
                    for i in range(1, 23): 
                        if selected_metrics == str(i):
                            st.markdown("**T-WAS Control " + str(i) + " Result**")
                            st.table(twas_control[str(i)])
                    st.markdown("***")

                    #tewas
                    for j in range(1, 23):         
                        if selected_metrics_1 == str(j):
                            st.markdown("**TE-WAS Control " + str(j) + " Result**")
                            st.table(tewas_control[str(j)])  
                    st.markdown("***")

                    #pewas
                    for k in range(1, 23):         
                        if selected_metrics_2 == 'CSF':
                            if selected_metrics_3 == 'Cardiometabolic': 
                                if selected_metrics_4 == str(k):
                                    try: 
                                        st.markdown("**P-WAS CSF Cardiometabolic Control " + str(k) + " Result**")
                                        st.table(pwas_csf_cardio_control[str(k)]) 
                                    except:
                                        #continue
                                        st.table(df_1)
                        if selected_metrics_2 == 'CSF':
                            if selected_metrics_3 == 'Inflammation': 
                                if selected_metrics_4 == str(k):
                                    try:
                                        st.markdown("**P-WAS CSF Inflammation Control " + str(k) + " Result**")
                                        st.table(pwas_csf_inf_control[str(k)])
                                    except:
                                        st.table(df_1)
                        if selected_metrics_2 == 'CSF':
                            if selected_metrics_3 == 'Neurology': 
                                if selected_metrics_4 == str(k):
                                    try:
                                        st.markdown("**P-WAS CSF Neurology Control " + str(k) + " Result**")
                                        st.table(pwas_csf_neuro_control[str(k)])
                                    except:
                                        st.table(df_1)
                        if selected_metrics_2 == 'CSF':
                            if selected_metrics_3 == 'Oncology': 
                                if selected_metrics_4 == str(k):
                                    try:
                                        st.markdown("**P-WAS CSF Oncology Control " + str(k) + " Result**")
                                        st.table(pwas_csf_onc_control[str(k)])
                                    except:
                                        st.table(df_1)
                        if selected_metrics_2 == 'Plasma':
                            if selected_metrics_3 == 'Cardiometabolic': 
                                if selected_metrics_4 == str(k):
                                    try:
                                        st.markdown("**P-WAS Plasma Cardiometabolic Control " + str(k) + " Result**")
                                        st.table(pwas_plasma_cardio_control[str(k)]) 
                                    except:
                                        st.table(df_1)
                        if selected_metrics_2 == 'Plasma':
                            if selected_metrics_3 == 'Inflammation': 
                                if selected_metrics_4 == str(k):
                                    try:
                                        st.markdown("**P-WAS Plasma Inflammation Control " + str(k) + " Result**")
                                        st.table(pwas_plasma_inf_control[str(k)])  
                                    except:
                                        st.table(df_1)
                        if selected_metrics_2 == 'Plasma':
                            if selected_metrics_3 == 'Neurology': 
                                if selected_metrics_4 == str(k):
                                    try:
                                        st.markdown("**P-WAS Plasma Neurology Control " + str(k) + " Result**")
                                        st.table(pwas_plasma_neuro_control[str(k)])
                                    except:
                                        st.table(df_1) 
                        if selected_metrics_2 == 'Plasma':
                            if selected_metrics_3 == 'Oncology': 
                                if selected_metrics_4 == str(k):
                                    try:
                                        st.markdown("**P-WAS Plasma Oncology Control " + str(k) + " Result**")
                                        st.table(pwas_plasma_onc_control[str(k)])  
                                    except:
                                        st.table(df_1)

                ########################  left column   #########################################    

                with left_column:
                #    fig = go.Figure()
                    #twas
                    for i in range(1, 23): 
                        if selected_metrics == str(i):
                            st.markdown("**T-WAS Case " + str(i) + " Result**")
                            st.table(twas_case[str(i)])   
                    #st.text("")
                    st.markdown("***")

                    #tewas
                    for j in range(1, 23): 
                        if selected_metrics_1 == str(j):
                            st.markdown("**TE-WAS Case " + str(j) + " Result**")
                            st.table(tewas_case[str(j)]) 
                    st.markdown("***")  

                    #pwas
                    for k in range(1, 23):
                        if selected_metrics_2 == 'CSF':
                            if selected_metrics_3 == 'Cardiometabolic':
                                if selected_metrics_4 == str(k):
                                    try:
                                        st.markdown("**P-WAS CSF Cardiometabolic Case " + str(k) + " Result**")
                                        st.table(pwas_csf_cardio_case[str(k)]) 
                                    except:
                                        #continue
                                        st.table(df_1)
                        if selected_metrics_2 == 'CSF':
                            if selected_metrics_3 == 'Inflammation': 
                                if selected_metrics_4 == str(k):
                                    try:
                                        st.markdown("**P-WAS CSF Inflammation Case " + str(k) + " Result**")
                                        st.table(pwas_csf_inf_case[str(k)]) 
                                    except:
                                        st.table(df_1)
                        if selected_metrics_2 == 'CSF':
                            if selected_metrics_3 == 'Neurology': 
                                if selected_metrics_4 == str(k):
                                    try:
                                        st.markdown("**P-WAS CSF Neurology Case " + str(k) + " Result**")
                                        st.table(pwas_csf_neuro_case[str(k)]) 
                                    except:
                                        st.table(df_1)
                        if selected_metrics_2 == 'CSF':
                            if selected_metrics_3 == 'Oncology': 
                                if selected_metrics_4 == str(k):
                                    try: 
                                        st.markdown("**P-WAS CSF Oncology Case " + str(k) + " Result**")
                                        st.table(pwas_csf_onc_case[str(k)]) 
                                    except:
                                        st.table(df_1)
                        if selected_metrics_2 == 'Plasma':
                            if selected_metrics_3 == 'Cardiometabolic': 
                                if selected_metrics_4 == str(k):
                                    try:
                                        st.markdown("**P-WAS Plasma Cardiometabolic Case " + str(k) + " Result**")
                                        st.table(pwas_plasma_cardio_case[str(k)]) 
                                    except:
                                        st.table(df_1)
                        if selected_metrics_2 == 'Plasma':
                            if selected_metrics_3 == 'Inflammation': 
                                if selected_metrics_4 == str(k):
                                    try:
                                        st.markdown("**P-WAS Plasma Inflammation Case " + str(k) + " Result**")
                                        st.table(pwas_plasma_inf_case[str(k)])  
                                    except:
                                        st.table(df_1)
                        if selected_metrics_2 == 'Plasma':
                            if selected_metrics_3 == 'Neurology': 
                                if selected_metrics_4 == str(k):
                                    try:
                                        st.markdown("**P-WAS Plasma Neurology Case " + str(k) + " Result**")
                                        st.table(pwas_plasma_neuro_case[str(k)]) 
                                    except:
                                        st.table(df_1) 
                        if selected_metrics_2 == 'Plasma':
                            if selected_metrics_3 == 'Oncology': 
                                if selected_metrics_4 == str(k):
                                    try:
                                        st.markdown("**P-WAS Plasma Oncology Case " + str(k) + " Result**")
                                        st.table(pwas_plasma_onc_case[str(k)])  
                                    except:
                                        st.table(df_1)



                # The plot!

                ternary_fig = px.scatter_ternary(df, a="T-WAS", b="P-WAS", c="TE-WAS", hover_name="Gene Symbol | Ensemble", 
                                                 color="Mean Z", color_continuous_scale=px.colors.sequential.Agsunset, opacity=0.5 )
                ternary_fig.update_traces(mode='markers', marker_line_width=2, marker_size=10)
                ternary_fig.update_layout(width=1400, height=800)
                #ternary_fig.update_layout(title='Ternary plot of potential multi-omic connections')

                st.markdown("**Ternary plot of potential multi-omic connections**")
                st.write(ternary_fig)

                ###Add ternary plot for case vs control
                #Dan's idea
                ###might be cool to have small manhattan plots per gene, showing eqtl locations with pvals on the multiomics app

              
            else:
                st.warning("Incorrect Username/Password")
            
    elif choice == "SignUp":
        st.subheader("Create New Account")
        new_user = st.text_input("Username")
        new_password = st.text_input("Password", type='password')
        
        if st.button("Signup"):
            create_usertable()
            add_userdata(new_user, new_password)
            st.success("You have successfully created a valid Account")
            st.info("Go to Login Menu to login")

        
        
        
        
        
        
        
        
if __name__ == '__main__':
    main()
