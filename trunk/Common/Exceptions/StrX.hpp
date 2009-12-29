/*************************************************************************************
    This is a simple class that lets us do easy (though not terribly efficient)
    trancoding of XMLCh data to local code page for display to report errors
    in cerr << line

	Dependencies
		none
*************************************************************************************/

#ifndef _STRX_

#define _STRX_

class StrX
{
    public :
        //  Constructor
        StrX(const XMLCh* const toTranscode)
        {
            // Call the private transcoding method
            fLocalForm = XMLString::transcode(toTranscode);
        }

        // Destructor
        ~StrX()
        {delete [] fLocalForm;}

        //  Getter method
        const char* localForm() const
        {return fLocalForm;}
    
    private :
        char*   fLocalForm;
};

inline ostream& operator<<(ostream& target, const StrX& toDump)
{
    target << toDump.localForm();
    return target;
}

#endif
