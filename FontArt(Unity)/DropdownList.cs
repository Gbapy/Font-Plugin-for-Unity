using System;
using System.Collections.Generic;
using UnityEngine;
#if UNITY_EDITOR
using UnityEditor;
#endif

public class DropdownList : PropertyAttribute
{
    public delegate string[] GetStringList();
    public DropdownList(string[] list)
    {
        List = list;
    }

    public DropdownList(Type type, string methodName)
    {
        var method = type.GetMethod(methodName);
        if (method != null)
        {
            List = method.Invoke(null, null) as string[];
        }
        else
        {
            Debug.LogError("NO SUCH METHOD " + methodName + " FOR " + type);
        }
    }
    public string[] List
    {
        get;
        private set;
    }
}
