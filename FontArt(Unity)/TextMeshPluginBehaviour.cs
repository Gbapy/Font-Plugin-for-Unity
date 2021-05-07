using System;
using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

[ExecuteInEditMode]
public class TextMeshPluginBehaviour : MonoBehaviour
{
#if UNITY_EDITOR
    //Plugin function that generates the 3D text model according to parameters in inspector.
    [DllImport("libFontArt", EntryPoint = "?generateModel@@YAHPEADH0H_N11MNNPEAH@Z")]
    private static extern int generateModel(UIntPtr text, int textSize, UIntPtr fontName, int fontNameSize
        , bool bBold, bool bUnderline, bool bItalic, float thickness, double xTwist, double yTwist, UIntPtr bufferInfo);

    //Plugin function that retrieves the 3D text mesh data.
    [DllImport("libFontArt", EntryPoint = "?getMeshData@@YAXPEAH0PEAM01@Z")]
    private static extern void getMeshData(UIntPtr vBufInfo, UIntPtr mBufInfo, UIntPtr vBuf, UIntPtr mBuf, UIntPtr nBuf);



    [DropdownList(typeof(PropertyDrawerHelper), "AllSystemFonts")]
    public string font = "Broadway";                                //Font of 3D Text model
    [TextArea]
    public string text = "FontArt For Unity3D!";                    //Content of 3D Text model
    public bool bold;                                               //Bold font
    public bool italic;                                             //Italic font
    public float thickness = 1;                                     //Thickness of 3D Text model
    public Vector2 curving;                                         //curving of 3D Text model

    private bool bParamChanged = false;                             //Flag according to changes in inspector
    // Start is called before the first frame update
    void Start()
    {
        if(GetComponent<MeshFilter>() == null)
        {
            gameObject.AddComponent<MeshFilter>();      //Add MeshFilter component.
        }
        if(GetComponent<MeshRenderer>() == null)
        {
            gameObject.AddComponent<MeshRenderer>();    //Add MeshRenderer component.
        }
    }

    private void OnValidate()
    {
        bParamChanged = true;                           //Flag according to changes in inspector
    }
    // Update is called once per frame
    void Update()
    {
        if (bParamChanged == false) return;
        unsafe
        {
            int ret = 0;                                                    //Return Value.
            byte[] txtBuffer = System.Text.Encoding.ASCII.GetBytes(text);   //Byte array from the string.
            byte[] fntBuffer = System.Text.Encoding.ASCII.GetBytes(font);   //Byte array from the font name.
            int[] bufferInfo = new int[3];                                  //Temperary buffer (used as return value).

            fixed(byte *txtBuf = txtBuffer, fntBuf = fntBuffer)
            {
                fixed(int *bufInfo = bufferInfo)
                {
                    ret = generateModel(new UIntPtr(txtBuf), text.Length, new UIntPtr(fntBuf), font.Length, bold, false, italic, 
                        thickness, curving.x, curving.y, new UIntPtr(bufInfo)); //Generate a 3D model

                }
            }
            if (ret != 0)
            {                                                       //Error Happened!
                Debug.LogWarning("Error happened!");                //Log Warning.
                goto _final;                                        //Final operation.
            }

            int[] vBufferInfo = new int[bufferInfo[0]];             //Integer array for vertex count.
            int[] mBufferInfo = new int[bufferInfo[0]];             //Integer array for triangle count.
            float[] vBuffer = new float[bufferInfo[1] * 3];         //Float array for the vertex.
            float[] nBuffer = new float[bufferInfo[1] * 3];         //Float array for the normals.
            int[] mBuffer = new int[bufferInfo[2] * 3];             //Integer array for the triangles.

            fixed(int *vBufInfo = vBufferInfo, mBufInfo = mBufferInfo, mBuf = mBuffer)
            {
                fixed(float *vBuf = vBuffer, nBuf = nBuffer)
                {
                    getMeshData(new UIntPtr(vBufInfo), new UIntPtr(mBufInfo), new UIntPtr(vBuf), 
                        new UIntPtr(mBuf), new UIntPtr(nBuf));      //Retrieve mesh data from plugin.
                }
            }

            Mesh mesh = new Mesh();                                 //Create a Mesh instance.
            List<Vector3> verts = new List<Vector3>();              //Array for vertex.
            List<Vector3> norms = new List<Vector3>();              //Arrat for vertex normal.
            List<int> tries = new List<int>();                      //Array for triangles.
            int vCount = 0;
            int mCount = 0;
            int p = 0;
            for (int i = 0; i < bufferInfo[0]; i++) {
                for(int j = 0;j < vBufferInfo[i];j++)
                {
                    verts.Add(new Vector3(vBuffer[vCount * 3], vBuffer[vCount * 3 + 1], vBuffer[vCount * 3 + 2]));
                    norms.Add(new Vector3(nBuffer[vCount * 3], nBuffer[vCount * 3 + 1], nBuffer[vCount * 3 + 2]));
                    vCount++;
                }
            }
            mesh.SetVertices(verts);                                //Set vertex.
            mesh.SetNormals(norms);                                 //Set normals
            
            for (int i = 0; i < bufferInfo[0]; i++)
            {
                //int[] tries = new int[mBufferInfo[i] * 3];
                for (int j = 0; j < mBufferInfo[i]; j++)
                {
                    tries.Add(p + mBuffer[mCount * 3]);
                    tries.Add(p + mBuffer[mCount * 3 + 1]);
                    tries.Add(p + mBuffer[mCount * 3 + 2]);
                    mCount++;
                }
                p += vBufferInfo[i];
            }
            mesh.SetTriangles(tries, 0);                            //Set triangles

            //Clear previous mesh.
            if (GetComponent<MeshFilter>().sharedMesh != null) GetComponent<MeshFilter>().sharedMesh.Clear();
            GetComponent<MeshFilter>().sharedMesh = mesh;           //Set Mesh.

            if (GetComponent<MeshRenderer>().sharedMaterial == null)
            {
                Material mat = new Material(Shader.Find("Standard"));
                mat.color = new Color(0, 0.5f, 0);
                GetComponent<MeshRenderer>().sharedMaterial = mat;  //Set material
            }
        }
    _final:
        bParamChanged = false;
    }
#endif
}
