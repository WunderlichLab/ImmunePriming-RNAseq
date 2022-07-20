import pandas as pd
import time
import os
import sys
from pathlib import Path


def get_ijwd():
    """Retrieve wd organizer.py is in, return path as Path object."""

    path = Path(os.getcwd())
    return path


def list_txt(dir_w_txts):
    """Take Path object with ImageJ txt files, return list only txt files in Path() form."""

    txtlist = []
    # Iterate through directory
    for o in dir_w_txts.iterdir():
        
        # only add .txt files
        if os.path.splitext(str(o))[-1].lower() == '.txt':
            txtlist.append(o)

    return txtlist


def list_dir(dir_total):
    """Retrieve all files in current directory."""

    dir_list = []

    for o in dir_total.iterdir():
        dir_list.append(o)

    return dir_list


def init_dataframe():
    """Create empty dataframe."""


    data = {"Genotype": [],
               "Condition": [],
               "TimePoint": [],
               "InjectionDate": [],
               "Replicate": [],
               "DF": [],
               "ColonyCount": [],
               "CFU_DF": [],
               "CFU_Fly": []}

    df = pd.DataFrame(data, columns = ["Genotype", "Condition", "TimePoint", "InjectionDate", "Replicate", "DF", "ColonyCount", "CFU_DF", "CFU_Fly"])

    return df


def init_cfuframe():
    """Create empty dataframe."""

    data = {"Genotype": [],
               "Pathogen": [],
               "Conc": [],
               "Hours": [],
               "Days": [],
               "CFU": [],
               "InjDate": [],
               "Notes": []}

    df = pd.DataFrame(data, columns = ["Genotype", "Pathogen", "Conc", "Hours", "Days", "CFU", "InjDate", "Notes"])

    return df


def fill_cfu(frame, infolist):
    """Fill information for cfu.csv."""
    
    genotype = [infolist[0][0]]
    pathogen = ["Efae"]
    conc = [infolist[1][0]]
    hours = [0]
    days = [infolist[2][0]]
    CFU = [infolist[7][0]]
    injdate = [infolist[3][0]]
    notes = ['']

    append_df = pd.DataFrame({"Genotype": genotype,
               "Pathogen": pathogen,
               "Conc": conc,
               "Hours": hours,
               "Days": days,
               "CFU": CFU,
               "InjDate": injdate,
               "Notes": notes})
    fill_df = frame.append(append_df)

    return fill_df


def fill_f(empty_frame, plst, w_dir, flist):
    """Populate the empty dataframe with file information and prompt for extra information."""

    txtlist = []
    replist = []
    DFlist = []
    CClist = []

    #print(flist)
    if flist == '':
            
        # Retrieve Information
        df_number = input("What factor are the dilutions? Enter an integer: ")

        # Check that dilution factor entered is an integer
        try:
            df_number = int(df_number)
        except:
            print("ERROR, not an integer")

        while type(df_number) != int:
            df_number = input("Enter an integer: ")
            try:
                df_number = int(df_number)
            except:
                print("ERROR, not an integer")
        genotype = str(input("What is the genotype?: "))
        condition = str(input("What is the condition?: "))
        assaydate = str(input("When is the time point?('d_'): "))
        injectiondate = str(input("When was the injection date? (MM/DD/YY): "))
    else:
        df_number = flist[6]

    # Retrieve specific txt file name as string, append to list
    for pth in plst:
        txt_string = str(os.path.split(pth)[-1])
        txtlist.append(txt_string)

    # Easier read paths for whole directory
    w_list = []
    for wpth in w_dir:
        wtxt_string = str(os.path.split(wpth)[-1])
        w_list.append(wtxt_string)
        
    # Check if .tif files have a .tif_results.txt file
    tif_lst = []
    
    for file in w_list:
        if os.path.splitext(str(file))[-1].lower() == '.tif':
            tif_lst.append(file)

    # Retrieve replicate number and calculate true DF
    for index, o in enumerate(tif_lst):
        tif_lst[index] = str(o)[0:3]
        rep = o[2]
        DF = pow(int(df_number), (int(o[0]) - 1))
        replist.append(rep)
        DFlist.append(DF)

    temp_lst = []
    for f in txtlist:
        temp_lst.append(str(f)[0:3])
    cut_lst = tif_lst.copy()
    for t in temp_lst:
        if t in cut_lst:
            cut_lst.remove(t)

    # Count number of colonies, request for zero or uncountable counts
    col_length = len(tif_lst)
    CClist = []
    
    for o in tif_lst:
        CClist.append(o)

    plist = []

    if len(temp_lst) != col_length:
        print("Missing .tif_results.txt file(s), answer following prompts")

    for index, c in enumerate(CClist):
        if c in cut_lst:
            U_Z = ''
            while U_Z != 'U' and U_Z != 'u' and U_Z != '0':
                U_Z = str(input("Is " + str(c) + " uncountable or zero?(U or 0): "))
            if U_Z == '0':
                U_Z = int(U_Z)
            if U_Z == 'u':
                U_Z = 'U'
            CClist[index] = U_Z

    for ind, p in enumerate(plst):
        i = str(os.path.split(p)[-1])
        i = i[0:3]
        if i in tif_lst:
            mod_ind = tif_lst.index(i)
            with open(p) as f:
                lines = f.readlines()
                CC = len(lines) - 1
                CClist[mod_ind] = CC

    
    # Calculate CFU per dilution and CFU per fly
    CFU_df_avg_list, CFU_fly_avg_list = cfu_calc(col_length, DFlist, CClist, df_number)
    
    # Repeating Lists for Genotype, Condition, Assay Date, and Injection Date for length matching
    glist  = []
    conlist = []
    adlist = []
    idlist = []

    c = 0
    while c < col_length and flist == '':
        glist.append(genotype)
        conlist.append(condition)
        adlist.append(assaydate)
        idlist.append(injectiondate)
        c += 1

    # Two flies, same metadata
    if flist != '':
        glist = flist[0]
        conlist = flist[1]
        adlist = flist[2]
        idlist = flist[3]
        replist = flist[4]
        DFlist = flist[5]
    
    append_df = pd.DataFrame({"Genotype": glist,
                            "Condition": conlist,
                            "TimePoint": adlist,
                            "InjectionDate": idlist,
                            "Replicate": replist,
                            "DF": DFlist,
                            "ColonyCount": CClist,
                            "CFU_DF": CFU_df_avg_list,
                            "CFU_Fly": CFU_fly_avg_list})
    fill_df = empty_frame.append(append_df)
    flist = [glist, conlist, adlist, idlist, replist, DFlist, df_number, CFU_fly_avg_list]

    return fill_df, flist


def cfu_calc(list_length, list_df, list_cc, df_num):
    """Output CFU per dilution and CFU per fly."""

    cfu_ind_list = []
    cfu_df_list = []
    cfu_fly_list = []

    # Populate empty lists
    c = 0
    while c < list_length:
        cfu_ind_list.append('')
        cfu_df_list.append('')
        cfu_fly_list.append('')
        c += 1

    
    # Calculate indiviudal CFU
    for index, indf in enumerate(list_df):
        cc = list_cc[index]

        if cc != 0 and cc != 'U' and cc != 'u':
            cfu = (cc * (indf / df_num)) * 250
            cfu_ind_list[index] = cfu
        elif cc == 'U':
            cfu_ind_list[index] = 'U'

    # print("CFU_ind_list: ", cfu_ind_list)


    # Calculate df group CFU
    for index, incfu in enumerate(cfu_ind_list):
        if index == 0 or index == 3 or index == 6 or index == 9 or index == 12 or index == 15 or index == 18 or index == 21:
            if index < list_length - 1:
                index_tuple = (cfu_ind_list[index], cfu_ind_list[index + 1], cfu_ind_list[index + 2])
                num_list = []
                for obj in index_tuple:
                    if type(obj) == float or type(obj) == int:
                        num_list.append(obj)
                    elif obj == 'U':
                        
                        num_list.append(obj)
        
                # print("NUM LIST: ", num_list)
                if 'U' in num_list:
                    nums = [i for i in num_list if type(i) != str]
                    if len(nums) != 0:
                        # Padding for calculating DF averages at higher dilutions, U is treated as missing datapoint, not used in average length
                        df_avg = sum(nums) / len(nums)
                    else:
                        df_avg = ""

                else:
                    # Padding for calculating DF averages at lower dilutions, 0 is treated as 0 datapoint, used in average length
                    df_avg = sum(num_list) / 3
        
                if len(num_list) == 0:
                    df_avg = ""
                cfu_df_list[index] = df_avg
                
    # print("num_list: ", num_list)
    # print("CFU_df_list: ", cfu_df_list)

    # Calculate fly CFU
    temp = []
    o = 0
    for index, dfcfu in enumerate(cfu_df_list):
        if type(dfcfu) == int or type(dfcfu) == float:
            num = float(dfcfu)
            temp.append(num)
            o += 1
    # print("temp: ", temp)
    avg = sum(temp) / o
    cfu_fly_list[0] = avg

    # print("CFU_fly_list: ", cfu_fly_list)
        
    return cfu_df_list, cfu_fly_list


def run():
    """Run program."""

    # TEST MODE
    test = 0

    # Two flies or one fly
    two = 0
    if test == 1:
        print("***TEST MODE***")
    two = str(input("Two Flies - 1 | One Fly - 0: "))
    while two != "0" and two != "1":
        two = str(input("Two Flies - 1 | One Fly - 0: "))
    two = int(two)
    print()
    flist = ''
    
    if two == 1:
        ij_txt_dir = get_ijwd()
        whole_dir = list_dir(ij_txt_dir)
        store_list = []
        cfu_list = []
        flist = ''
        dir_list = [m for m in whole_dir if m.is_dir()]
        i = 1
        
        for o in dir_list:
            print()
            print("Plate " + str(i))
            print()
            new = o
            new_d = list_dir(new)
            new_list = list_txt(new)
            new_blank = init_dataframe()
            new_frame, flist = fill_f(new_blank, new_list, new_d, flist)

            cfu_frame = init_cfuframe()
            cfu_fill = fill_cfu(cfu_frame, flist)

            store_list.append(new_frame)
            cfu_list.append(cfu_fill)
            i += 1

        #print(store_list)
        pop_frame = store_list[0].append(store_list[1])
        cfu_fill = cfu_list[0].append(cfu_list[1])
        
    else:
    # Get current directory with text files and create a data frame
        ij_txt_dir = get_ijwd()
        whole_dir = list_dir(ij_txt_dir)
        path_list = list_txt(ij_txt_dir)
        blank_frame = init_dataframe()
        pop_frame, flist = fill_f(blank_frame, path_list, whole_dir, flist)
        cfu_frame = init_cfuframe()
        cfu_fill = fill_cfu(cfu_frame, flist)

    # Ask for file name, default "" is Colonycounts.csv, file must end in .csv
    file_name = str(input("Enter existing or new .csv file name (Default: 'Colonycounts.csv'): "))
    if file_name == "":
        file_name = "Colonycounts.csv"

    while os.path.splitext(file_name)[-1].lower() != '.csv':
        file_name = str(input("Please enter file name with '.csv' extension at the end: "))
        if file_name == "":
            file_name = "Colonycounts.csv"
        
    csv_file = ij_txt_dir / file_name
    cfu_file = ij_txt_dir/ "cfu.csv"
    whole_dir = list_dir(ij_txt_dir)

    # Continue operation or quit program
    process = str(input("Continue - 1 | Quit - 0: "))
    while process != "1" and process != "0":
        process = str(input("Continue - 1 | Quit - 0: "))

    if process == '0':
        print("Quitting...")
        time.sleep(2)
        sys.exit(1)
    else:
        pass

    if test == 0:
        
    # Create file with specified name, or append to file if already in current working directory
        if csv_file not in whole_dir:
            print("Created " + file_name + "!")
            print("Created " + "cfu.csv" + "!")
            pop_frame.to_csv(csv_file, index = False)
            cfu_fill.to_csv(cfu_file, index = False)
        else:
            print("Appended data to " + file_name + "!")
            print("Appended data to " + "cfu.csv" + "!")
            pop_frame.to_csv(csv_file, mode = 'a', index = False, header = False)
            cfu_fill.to_csv(cfu_file, mode = 'a', index = False, header = False)

        print("Quitting...")
        time.sleep(2)
  
    else:
        print("***TEST MODE: Did not append or create file***")
        print("Quitting...")
        time.sleep(2)


if __name__ == '__main__':
    run()
