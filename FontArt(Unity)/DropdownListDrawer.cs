using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEditor;

#if UNITY_EDITOR
[CustomPropertyDrawer(typeof(DropdownList))]
public class StringInListDrawer : PropertyDrawer
{
    // Draw the property inside the given rect
    public override void OnGUI(Rect position, SerializedProperty property, GUIContent label)
    {
        DropdownList dpList = attribute as DropdownList;
        string[] list = dpList.List;
        if (property.propertyType == SerializedPropertyType.String)
        {
            int index = Mathf.Max(0, Array.IndexOf(list, property.stringValue));
            index = EditorGUI.Popup(position, property.displayName, index, list);
            property.stringValue = list[index];
        }
        else
        {
            base.OnGUI(position, property, label);
        }
    }
}
#endif