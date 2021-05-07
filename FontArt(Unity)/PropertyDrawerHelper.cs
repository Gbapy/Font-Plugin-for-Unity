using System;
using System.Runtime.InteropServices;
using System.Collections.Generic;
using UnityEditor;
using UnityEngine;

public static class PropertyDrawerHelper
{
#if UNITY_EDITOR
    //Plugin function for initialization.
    [DllImport("libFontArt", EntryPoint = "?initOperation@@YAHXZ")]
    private static extern int initOperation();

    //Plugin function that retrieves the list of font names.
    [DllImport("libFontArt", EntryPoint = "?getFontList@@YAHPEAD@Z")]
    private static extern int getFontList(UIntPtr buffer);

    public static string[] AllSystemFonts()
    {
        List<string> temp = new List<string>();         //String array for return.
        int fCount = initOperation();                   //Initialization of plugin, Getting the numbers of fonts installed in OS.
        if (fCount > 0)
        {
            byte[,] buffer = new byte[fCount, 255];     //Byte array for font names
            unsafe
            {
                fixed (byte* buf = buffer)
                {
                    getFontList(new UIntPtr(buf));      //Getting the raw font list
                }
            }
            for (int i = 0; i < fCount; i++)
            {
                byte[] b = new byte[buffer[i, 0]];      //Byte array for real font name      
                for (int j = 1; j <= buffer[i, 0]; j++)
                {
                    b[j - 1] = buffer[i, j];
                }
                string fName = System.Text.Encoding.ASCII.GetString(b);
                temp.Add(fName);                        //Add a real font name.
            }
            temp.Sort();                                //Sorting.
            return temp.ToArray();                      
        }
        else
        {
            Debug.LogWarning("Unknown error happened when getting system font list");
        }
        return null;
    }

#endif
}