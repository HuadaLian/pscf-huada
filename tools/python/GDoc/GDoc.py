#!/usr/bin/env python
#****M root/GDoc ------------------------------------------------------
# MODULE
#    GDoc - Generic documentation tool
# PURPOSE
#    GDoc is designed to extract multi-line documentation strings from
#    programs or other text files, and output them formatted documentation.
#    Currently, the program can produce a set of hyperlinked html files, 
#    a single latex document, or (for storage) xml documents.
# FORMAT 
#    GDoc extracts documentation "Items", which are usually header blocks
#    written above functions or variables, and outputs these as formatted 
#    documentation. Each documentation item contains a series of sections,
#    each of which contains a heading on a line by itself (e.g., "PURPOSE",
#    "ARGUMENTS", etc.) followed by one or more lines of text. To provide
#    an example, the source code for the GDoc module is itself formatted 
#    with header blocks in the required formate.
#
#    Items: 
#        GDoc recognizes the beginnning and end of each documentation 
#    Item by using the Python regular expression engine to search for 
#    beginning and ending marker lines. Examples of such marker lines
#    are show below:
#
#    beginning marker:    "#****c  parent/name" , or
#                         "#****ic parent/name"
#
#    ending marker        "#***"
#   
#    Here # is an example of a comment character that is used to mark the 
#    beginning of a one line comment or remark in the language used in the 
#    source code (e.g., ! in Fortran 90). The 4 asterisks in the beginning 
#    marker line and the 3 asterisks in the ending marker are literal and 
#    are required.  The "c" following the asterisks in the first beginning
#    marker is an example of a single character code used to indicate the 
#    type of object being documented (e.g., "c" for class, "f" for function,
#    etc.) in a user defined dictionary. The optional "i" character in the 
#    second example of a beginning marker is used to indicate that the 
#    item is an "internal" (or private) item.  The string "name" in the above 
#    denotes the name of the documentation Item, which is typically the name 
#    of the variable, function, or other thing that is being documented.  
#    GDoc searches for instances of Item names to create hyperlinks in html 
#    documentation.  The string "parent" in the above is the name of 
#    another documentation Item that is the parent of this one in a hierarchy: 
#    For example, a module may be the parent of a variable or procedure, and 
#    a class may be the parent of its data attributes and/or methods.  These 
#    parent/child relationships are used to create a hierarchical tables of 
#    contents for each source file (in multi-file html documentation) and 
#    for the entire document (in both html or latex documentation). 
#
#    Sections:
#         Each documentation item contains a list of Sections, each
#    of which begins with a heading, such "NAME" or "PURPOSE", on a line
#    by itself, followed by one or more lines of text. The user must
#    provide supply a list of recognized Section headings in the 
#    configuration file, as discussed below.
#
#    Configuration:
#         Before parsing the source files, GDoc reads a configuration 
#    file. The configuration file must specify, among other things: 
#    i)  The comment character that appears at the beginning of 
#        begining and ending marker lines, and at the beginning of each
#        line of documentation (e.g, '#' in the above examples and in
#        the source code of the GDoc module).
#    ii) A list of names of the types of items to document (e.g., 
#        "module", "function", "variable", etc.), and an associated 
#        single character codes for each. These (key : value) pairs 
#        are stored in a user defined dictionary, e.g., 
#        ("m":"module", #    "f":"function", "v":"variable" , ... ).  
#    iii)A list of recognized section headings (e.g., "NAME", "PURPOSE", 
#        etc.).  
#   The format of the configuration file is discussed in detail in the 
#   documentation for function read_config. Examples of working 
#   configuration files for generating the hmtl and latex documentation 
#   for GDoc itself are given in the html.rc and latex.rc files, 
#   respectively.
# USAGE
#    To execute GDoc in Unix or Linux:
#
#    1) Change directory to a directory one level above the directory
#    containing the source file source, and the directory in which you
#    intend to create documentation.  Make the documentation directory,
#    if it does not already exist, and create a configuration file. 
#
#    2) Start the python interpreter by typing "python" at the shell
#    prompt. Within the interpreter, type "import GDoc". Then, enter 
#    the name of the configuration file (relative to the current 
#    directory) when prompted to do so. The GDoc directory must either 
#    be in the PYTHONPATH environment variable or a subdirectory of a 
#    directory in PYTHONPATH in order for python to find the module.
#   
#    To execute GDoc in MS Windows, one may either:
#
#    1) Start up the python interpreter and "import GDoc", as in
#    unix. This requires that python search path be set so that it 
#    can find the GDoc directory.
#
#    2) Double click on the GDoc/GDoc.py file, which should look 
#    like a snake. This, however, makes the GDoc directory the 
#    current directory, and thus requires that the source and 
#    documentation directory be specified in the configuration file
#    either as as absolute paths or as paths relative to the GDoc 
#    directory. 
#
#*** -------------------------------------------------------------------
import string, re
import os, shutil
import file_util

# Defaults: 
rc_filename = "GDoc.rc"  # default configuration file name
com_char    = "#"        # default comment character
    
# String containing characters allowed in variable names
name_char = string.letters + string.digits + '_' + '\.'

#****c GDoc/Section -----------------------------------------------------
# CLASS
#     Section
# PURPOSE
#     Contains one Section within an Item of documentation. Each
#     Section contains a heading (e.g., 'FUNCTION', 'PURPOSE', etc.)
#     and the corresponding text.
# CLASS ATTRIBUTES
#     types  = list of recognized heading strings, e.g.,
#              Section.types = ('FUNCTION','PURPOSE',..)
#     ignore = list of Section types to exclude from documentation
# INSTANCE ATTRIBUTES
#     heading  = heading string (must be in list Section.types)
#     text     = list of strings, each one line with '\n' stripped
# METHODS
#     __init__([heading],[text]) - constructor
#     write(format,namespace)    - returns string representation 
#                                  in format = 'html' or 'xml'
#*** ----------------------------------------------------------------------
class Section:

    types = []  # List of allowed Section heading strings, 
                # e.g., 'PURPOSE', 'ARGUMENTS', 'AUTHOR', etc.

    ignore = [] # List of heading strings for sections that should
                # NOT be included in documentation

    def __init__(self, heading='',text=[] ):
        self.heading = heading
        self.text    = text
    #end method Section.__init__


    #****m Section/Section.write ---------------------------------------
    # METHOD
    #     Section.write(format,namespace)
    # PURPOSE
    #     Output one documentation Section as a string
    # ARGUMENTS
    #     format    - output format ('html' or 'xml')
    #     namespace - dictionary of Items accessed by name. 
    #                 Used to hyperlink html output
    # RETURN
    #     string representation of item in chosen output format
    #     If format == 'html', result is hyper-linked using namespace
    #     If Section.type is in Section.ignore, returns blank string
    #*** --------------------------------------------------------------
    def write(self,format='xml',namespace=''):
        if self.heading in Section.ignore:
            return ''
        else:
            s = []
            # Heading
            if format == 'html' :
                s.append( r"<p><strong>" )
                s.append( self.heading )
                s.append( r"</strong></p>  <pre>" + '\n' )
            elif format == 'tex' :
                s.append( self.heading )
                s.append( r'\begin{verbatim}' + '\n' )
            elif format == 'xml' :
                s.append( r'<Section name="'+self.heading+ r'">'+'\n')
            # Text
            for line in self.text:
                if format == 'html':
                    s.append(link_namespace(line,namespace) + '\n' )
                elif format == 'tex' :
                    s.append(line + '\n' )
                elif format == 'xml' :
                    s.append( line + '\n' )
            # Closing
            if format == 'html' :
                s.append( r"</pre>"  + '\n' )
            elif format == 'tex' :
                s.append( r'\end{verbatim}' + '\n' )
            elif format == 'xml' :
                s.append( r'</Section>' + '\n' )
            return string.join(s,'')
    #end method Section.write

#end class Section


#****c GDoc/Item ------------------------------------------------------------
# CLASS
#     Item
# PURPOSE
#     Each Item object contains the documentation for one
#     module, variable, function, class, etc.
# CLASS ATTRIBUTES
#     types   = Dictionary of allowed Item types, 
#               keys are user defined one-character strings
#               values are strings, e.g., FUNCTION, VARIABLE, etc.
#     current = Item currently being processed
# INSTANCE ATTRIBUTES
#     name    = name (string)
#     type    = type of object (string)
#     parent  = parent (Item object)
#     access  = 'public' or 'private' (string)
#     SrcFile = file from which it was extracted (SrcFile object)
#     pattern = regular expression object for name
#     Sections= list of Section objects
# METHODS
#     __init__([attributes])  - constructor, all attributes optional
#     read(lines,i,n)         - extract documentation for Item from source
#     write(format,namespace) - returns string representation in 'format' 
#     html_ref(ext,target)    - returns hyper-link to self
#     html_link(line,namespace) - returns line with self.names replaced
#                                 by href links
#*** --------------------------------------------------------------------------
class Item:

    types   = {}    
    current = None  

    def __init__( self, name='', type='', parent='', access='public', \
                  SrcFile = '',pattern=''):
        self.type    = type      # character flag for type
        self.parent  = parent    # parent of Item
        self.name    = name      # name of Item
        self.pattern = pattern   # regular expression for name
        self.access  = access    # 'public' or 'private'
        self.SrcFile = SrcFile   # enclosing source file Object 
        self.Sections = []       # list of Section objects
        self.children = []       # list of child Items
    #end method Item.__init__


    #****m Item/Item.read ----------------------------------------------
    # METHOD
    #     Item.read(lines,i,n)
    # PURPOSE
    #     Read one documentation Item object from a list of text lines
    # ARGUMENTS
    #     lines - list of text line strings
    #     i     - # of current line
    #     n     - total # of strings in lines
    # RETURN
    #     i     - current line, from which to begin reading next Item
    #*** --------------------------------------------------------------
    def read(self,lines,i,n):
        global name_char
        Item.current = self
        while 1:
            if not i < n :
                return -1
            match = Item.bgn.match(lines[i])
            if match:
    
                # Parse beginning line marker
                line1       = string.strip(match.group(0))
                self.type   = line1[-1]
                line2       = lines[i][match.end():]
                line2       = re.sub(r'^(\*)*\s*','',line2)
                match2      = re.match(r"^.*\/",line2)
                self.parent = string.strip(match2.group(0)[:-1])
                line3       = line2[match2.end():]
                match3      = re.match('['+name_char+r']*',line3)
                line3       = line3[:match3.end()]
                self.name   = string.strip(line3)
                if line1[-2] == 'i' :
                    self.access = 'private'
                else :
                    self.access = 'public'
                self.pattern= re.compile(self.name)
                self.Sections= []
                a_Section    = None
                
                s=[]
                s.append( r'</pre>')
                s.append( r'<h2 id="' + self.name + r'" >' )
                s.append( self.name )
                s.append( r'</h2> <pre>'  )
                lines[i] = s
                
                i = i + 1
                break
            else :
                lines[i] = escape(lines[i],'html')
                i = i + 1
            
        while 1:
            if not i < n :
                raise "Error: Reached EOF searching for end marker"
            match = Item.end.match(lines[i]) 
            if match :
                if a_Section :
                    self.Sections.append(a_Section)
                lines[i] = escape(lines[i],'html')
                i  = i + 1
                return i 
            else:
                line1 = Item.com.sub(' ',lines[i])
                line2 = string.strip(line1)
                if line2 in Section.types :
                    it = i
                    if a_Section :
                        self.Sections.append(a_Section)
                    a_Section = Section(heading=line2,text=[])
                elif a_Section:
                    if a_Section.heading == 'SOURCE':
                         if not i == it+1 :
                            a_Section.text.append(lines[i])
                    else:
                         a_Section.text.append(line1)
                lines[i] = escape(lines[i],'html')
                i = i + 1
        else:
            raise "Error: Unexpected completion of while loop "
    #end method Item.read


    #****m Item/Item.write ---------------------------------------------
    # METHOD
    #     Item.write(format,namespace,Doc_access,recursive)
    # PURPOSE
    #     Output one documentation item as a string
    # ARGUMENTS
    #     format    - output format ('html', 'tex', or 'xml')
    #     namespace - dictionary of Items accessed by name. 
    #                 Used to hyperlink html output
    #     Doc_access- 'public'  -> write only if public Items,
    #                 'private' -> write public or private
    # RETURN
    #     string representation of item in chosen output format
    #     If format == 'html', result is hyper-linked using namespace
    #*** --------------------------------------------------------------
    def write(self,format='html',namespace='',Doc_access='public',recursive=''):
        Item.current = self
        if self.access == 'public' or Doc_access == 'private' :
            s = []
            
            # Opening
            if format == 'html' :
                s.append( r'<h2 id="' + self.name + r'" >' )
                s.append( self.html_ref('_',"source") )
                # s.append( self.parent + '/' + self.name  )
                s.append( r'</h2>' + '\n' )
            elif format == 'tex' :
                if self.depth == 1:
                    s.append( r'\newpage'  + '\n' )
                    s.append( r'\section{' )
                elif self.depth == 2:
                    s.append( r'\subsection{' )
                elif self.depth == 3:
                    s.append( r'\subsubsection{' )
                if self.depth > 0 :
                    Item_title = Item.types[self.type] + '  ' + self.name
                    s.append( escape(Item_title,'tex') )
                    s.append( r'}' + '\n' )
            elif format == 'xml' :
                s.append( r'<Item ')
                s.append( r'type="'   + self.type   + '" ')
                s.append( r'parent="' + self.parent + '" ' )
                s.append( r'name="'   + self.name   + '" ')
                s.append( r'SrcFile="'+ self.SrcFile.name + '" ')
                s.append( r'>'+'\n')
                
            # Sections
            for a_Section in self.Sections :
                s.append( a_Section.write(format,namespace) )

            # Closing
            if format == 'xml' :
                s.append( r'</Item>' + '\n' )

            if recursive :
                for child in self.children:
                    s.append(child.write(format,namespace,Doc_access,recursive))

            return string.join(s,'')
        else:
            return ''
    #end method Item.write



    #****m Item/Item.html_ref ------------------------------------------
    # METHOD
    #     Item.html_ref([ext],[target])
    # PURPOSE
    #     Returns string containing html element <a ... >Item.name</a>
    #     with an html href attribute linking to anchor of Item 
    # ARGUMENTS
    #     ext    = if present, appends 'ext' string to filename in href.
    #              Used to distinguish html documentation and html source
    #     target = if present, adds target attribute to href.
    #              Used to create new Browswer window
    #*** ---------------------------------------------------------------
    def html_ref(self,ext='',target=None):
        s = self.SrcFile.name + ext
        if SrcFile.current :
            if s == SrcFile.current.name :
                s = ''
            else:
                s = s + '.html'
                s = file_util.relative_path(SrcFile.current.name, s)
        else:
            s = s + '.html'
        s = '"' + s + '#' + self.name +'"'
        if target:
            s = s + ' target="' + target + '"'
        return r'<a href=' + s + r'>'+ self.name + r'</a>'

    
    #****m Item/Item.html_link -----------------------------------------
    # METHOD
    #     Item.html_link(line)
    # RETURN
    #     string line with instances of self.name replaced by html
    #     <a href=..>self.name</a> hyperlinks to self. Uses method html_ref
    #*** ---------------------------------------------------------------
    def html_link(self,line):
        name    = self.name
        pattern = self.pattern
        match   = re.search(pattern,line)
        if match:
            
            s = match.start()
            e = match.end()
            f = 1
            
            if s > 0:
                pre = line[:s]
                if line[s-1] in name_char : f = 0
            else:
                pre = ''
                
            if e < len(line) :
                post = line[e:]
                if line[e] in name_char : f = 0
            else:
                post = ''
                
            if f:
                if name == Item.current.name :
                    repl = r'<strong>' + name + r'</strong>'
                else:
                    repl =self.html_ref()
            else:
                repl = name
                
            line  = pre + repl + self.html_link(post)
        return line
    #end method link


    #****m Item/Item.set_depth -----------------------------------------
    # METHOD
    #     Item.set_depth(depth)
    # PURPOSE
    #     Sets depth attribute of self and all of its descendants.
    #***  --------------------------------------------------------------
    def set_depth(self,depth):
        self.depth = depth
        for child in self.children:
            child.set_depth(depth+1)

    
    #****m Item/Item.html_toc --------------------------------------------
    # METHOD
    #     Item.html_toc()
    # PURPOSE
    #     Returns string containing an html table of contents list item 
    #     for self, and a nested list containing its children, if any. 
    #     The method is applied recursively to children to generate a 
    #     table of contents for the tree of Items rooted at self. 
    #*** -----------------------------------------------------------------
    def html_toc(self):
        s = []
        if not self.type == 'root':
            s.append(r'<li>' + Item.types[self.type] + ' ')
            s.append(self.html_ref() + r'</li>' + '\n')
        if self.children:
            s.append( r'<ul>' + '\n' )
            for child in self.children:
                s.append(child.html_toc())
            s.append( r'</ul>' + '\n')
        return string.join(s,'')

#end Class Item


#****c GDoc/SrcFile ---------------------------------------------------
# CLASS
#     SrcFile - documentation extracted from a source file
# CLASS ATTRIBUTES
#     current = current SrcFile
# INSTANCE ATTRIBUTES
#     name  = SrcFile name (in format used for output file names)
#     Items = list of Item objects
#*** ------------------------------------------------------------------
class SrcFile:

    current = None

    def __init__( self, name='' ):
        self.name  = re.sub(r'\.','_',name) # SrcFile name
        self.Items = []                     # List of Items

    #****m SrcFile/SrcFile.read ---------------------------------------
    # METHOD
    #     SrcFile.read(file,namespace)
    # PURPOSE
    #     Extract documentation Items from an entire file
    # ARGUMENTS
    #     file      - file object
    #     namespace - namespace dictionary to which Items are added
    # RETURN
    #     lines  - Full text of input source file, with minimal
    #              hyper-text anchors inserted for documented items
    #*** -------------------------------------------------------------
    def read(self,file,namespace=''):
        SrcFile.current = self
        # Read File and chomp end-of-lines
        lines = file.readlines()  
        n = len(lines)
        for i in range(n):
            lines[i] = string.replace(lines[i],'\n','')  
            n = len(lines)
        i = 0
        while i < n :
            an_Item = Item(SrcFile=self)
            i = an_Item.read(lines,i,n)
            if not i == -1 :
                self.Items.append(an_Item)
                namespace[an_Item.name] = an_Item
            else:
                break
        return lines
    #end method SrcFile.read


    #****m SrcFile/SrcFile.write --------------------------------------
    # METHOD
    #     SrcFile.write(format,namespace,Doc_access)
    # PURPOSE
    #     Output documentation for instance SrcFile
    # ARGUMENTS
    #     format    - output format ('html' or 'xml')
    #     namespace - dictionary of Items accessed by name. 
    #                 Used to hyperlink html output
    #     Doc_access- 'public'  -> write only if public Items,
    #                 'private' -> write public or private
    # RETURN
    #     string representation of SrcFile in chosen output format
    #     If format == 'html', output is hyper-linked using namespace
    #*** --------------------------------------------------------------
    def write(self,format='xml',namespace='',Doc_access='public'):
        SrcFile.current = self
        s = []  # Character list, assembled by appendage
        # Top Matter
        s.append( header(self.name,format) )
        if format == 'html' :
            s.append( Document.current.index_links(self.name) )
            s.append( r'<h3 align="center">TABLE OF CONTENTS</h3>'+'\n')
            s.append( r'<ul>' )
            for an_Item in self.Items:
                if an_Item.access == 'public' or Doc_access == 'private':
                    if an_Item.parent == Document.current.root.name :
                        s.append( an_Item.html_toc() )
                    else:
                        parent = namespace[an_Item.parent]
                        if not parent in self.Items:
                            s.append( an_Item.html_toc() )
        elif format == 'xml' :
            s.append( r'<SrcFile name="' + self.name + '">' + '\n' )
        # Body: Output Items
        for an_Item in self.Items:
            s.append( an_Item.write(format,namespace,Doc_access) )
        # Bottom Matter
        if format == 'xml' :
            s.append( r'</SrcFile>' )
        s.append( footer(format) )
        return string.join(s,'')
    # end method SrcFile.write

# end Class SrcFile


#****c GDoc/Document --------------------------------------------------
# CLASS
#     Document - documentation for an entire program or project
# INSTANCE ATTRIBUTES
#     src_dir  = path name for source directory
#     doc_dir  = path name for documentation directory
#     filenames = list of file names of source files in src_dir
#     format   = format of documentation, 'html' or 'xml'
#     access   = 'public' (public Items only) or 'private' (private also)
#     SrcFiles = list of SrcFile objects, in same order as filenames
#     namespace = dictionary of Item objects accessed by name
#     indices  = dictionary of lists of items, used to create indices.
#                Keys are type name, or "master" for master index.
#     root     = string identifier for "parent" of top level Items 
# METHODS
#     read()  - read source files in self.filenames from self.src_dir
#     make_indices() - make master indices and 
#     make_tree() - finds parents of each Item, thus creating a tree
#     html_index(key,columns) - returns html table of self.indices[key]
#     index_links() - return string with hyperlinks to toc and indices
#     write() - write documentation to files in self.doc_dir
#*** ------------------------------------------------------------------
class Document:

    current = None

    def __init__(self,src_dir='',doc_dir='',filenames=[],\
                 format='xml',access='public'):
        self.src_dir   = src_dir
        self.doc_dir   = doc_dir
        self.filenames = filenames
        self.format    = format
        self.access    = access
        self.SrcFiles  = []
        self.namespace = {}
        self.indices   = {}
        self.root      = None

        
    #****m Document/Document.read --------------------------------------
    # METHOD
    #     Document.read()
    # PURPOSE
    #     Extract documentation for a multi-file Document
    #*** ---------------------------------------------------------------
    def read(self):
        Document.current = self
        os.chdir(self.src_dir)
        for filename in self.filenames:
        
            # Read a source file 
            print "Reading file " + filename
            file  = open(filename,'r')
            a_SrcFile = SrcFile(name=filename)
            source = a_SrcFile.read(file,self.namespace)
            self.SrcFiles.append(a_SrcFile)
        
            # If html format, echo source file as html in doc_dir
            if self.format == 'html' :
                htm_src = os.path.join(self.doc_dir,a_SrcFile.name+'_.html')
                file = file_util.open_w(htm_src)
                s=[]
                s.append(header(a_SrcFile.name,'html'))
                s.append(r'<pre>' + '\n')
                for line in source:
                    s.append( string.join(line,'') + '\n' )
                s.append(r'</pre>')
                file.write( string.join(s,'') )

            del a_SrcFile
    #end Document.read
   
        
    #****m Document/Document.make_indices ------------------------------
    # METHOD
    #     Document.make_indices()
    # PURPOSE
    #     Make dictionary of indices Document.indices 
    #     The keys of Document.indices are Item types, and "master".
    #     The values are sorted lists of Items of each type, and all Items.
    # USAGE
    #     Call after read method
    #*** ---------------------------------------------------------------
    def make_indices(self):
        names = {}
        names['master'] = []
        for a_SrcFile in self.SrcFiles:
            for an_Item in a_SrcFile.Items :
                if an_Item.access == 'public' or self.access == 'private' :
                    names['master'].append(an_Item.name)
                    type = Item.types[an_Item.type]
                    if not type in names.keys() :
                        names[type] = []
                    names[type].append(an_Item.name)
        for key in names.keys():
            names[key].sort()
            self.indices[key]=[]
            for name in names[key]:
                self.indices[key].append(self.namespace[name])


    #****m Document/Document.make_tree ---------------------------------
    # METHOD
    #     Document.make_tree()
    # PURPOSE
    #     Construct tree of parent and child items by constructing
    #     Item.children list for each Item in namespace
    # USAGE
    #     Call after read method
    #*** ---------------------------------------------------------------
    def make_tree(self):
        for a_SrcFile in self.SrcFiles:
            for an_Item in a_SrcFile.Items :
                if an_Item.access == 'public' or self.access == 'private' :
                    if an_Item.parent in self.namespace.keys():
                        parent = self.namespace[an_Item.parent]
                        parent.children.append(an_Item)
                    elif not self.root :
                        self.root = Item(name=an_Item.parent,type='root')
                        self.root.children.append(an_Item)
                    elif an_Item.parent == self.root.name :
                        self.root.children.append(an_Item)
                    else:
                        raise "Error in Document.make_tree" + '\n' \
			      "More than one root:" + '\n' \
			      + self.root.name + ' and ' + an_Item.parent
        # Set depth attributes for all Items in tree
        self.root.set_depth(0)
                
    #****m Document/Document.html_index --------------------------------
    # METHOD
    #     Document.html_index(key,columns)
    # PURPOSE
    #     Returns a string containing the html table containing an index
    # ARGUMENTS
    #     key       - key to index of interest in dictionary self.indices
    #     columns   - number of columns in html table
    #*** ---------------------------------------------------------------
    def html_index(self,key,columns):
        SrcFile.current = None
        index = self.indices[key]
        s = []  # Character list, assembled by appendage
        s.append( header(key,'html') )
        s.append( self.index_links() )
        s.append( r'<h1>' + key + ' index</h1>' + '\n' )
        s.append( r'<table cellspacing="6">' + '\n' )
        i = 0
        j = 0
        k = 0
        while i < len(index):
            an_Item = index[i]
            if j == 0:
                s.append( r'<tr>' + '\n' )
            s.append(r'<td class="index">'+an_Item.html_ref()+r'  </td>'+'\n')
            i = i + 1
            j = j + 1
            if j == columns :
                j = 0
                k = k + 1
                s.append( r'</tr>' + '\n' )
        if not j == 0 :
            s.append( r'</tr>' + '\n' )
        s.append( r"</table>" + '\n' )
        s.append( footer('html') )
        #s.append( r"</body>" + '\n' )
        #s.append( r"</html>"  )
        return string.join(s,'')

        
    #****m Document/Document.index_links -------------------------------
    # METHOD
    #     Document.index_links()
    # PURPOSE
    #     Returns a string containing hyper-links to indices and toc
    #*** ---------------------------------------------------------------
    def index_links(self, filename=''):
        s=[]
        s.append(r'<p>'+'\n')
        ref = file_util.relative_path(filename,"toc.html")
        ref = ref + "#top"
        s.append(r'[<a href="'+ref+r'">' + "table of contents" + r'</a>]'+'\n')
        ref = file_util.relative_path(filename,"master_index.html")
        ref = ref + "#top"
        s.append(r'[<a href="'+ref+r'">' + "master index" + r'</a>]'+'\n')
        for key in self.indices.keys():
            if not key == 'master' :
                button = string.strip(key)
                if button[-1] == 's' :
                    button = button + 'es'
                else :
                    button = button + 's'
                ref = file_util.relative_path(filename,key+"_index.html")
                ref = ref + "#top"
                s.append(r'[<a href="'+ref+r'">'+button+r'</a>]'+'\n')
        s.append(r'</p>'+'\n')
        return string.join(s,'')

    
    #****m Document/Document.write ------------------------------------
    # METHOD
    #     Document.write()
    # PURPOSE
    #     Output documentation for a multi-file Document
    # USAGE
    #     Call after read, make_indices, and make_tree methods
    #*** --------------------------------------------------------------
    def write(self):
        global GDoc_dir
        Document.current = self
        file_util.chdirs(self.doc_dir)

        if self.format == 'html' or self.format == 'xml' :

            for a_SrcFile in self.SrcFiles:
                filename = a_SrcFile.name
                file = file_util.open_w(filename + "." + self.format)
                file.write( a_SrcFile.write\
                            (self.format,self.namespace,self.access) )
                file.close()
                del a_SrcFile

        if self.format == 'html' :

            # Write index files
            for key in self.indices.keys():
                file = file_util.open_w(key + "_index.html")
                file.write(self.html_index(key,5))
                file.close()
    
            # Write toc file
            file = file_util.open_w("toc.html")
            file.write( header('Table of Contents','html') )
            file.write(self.index_links())
            file.write(r'<h3 align="center">TABLE OF CONTENTS</h3>'+'\n')
            file.write(self.root.html_toc())
            file.close()

        if self.format == 'tex':

            filename = raw_input\
               ("Enter filename (without .tex extension) for latex output: ")
            file = file_util.open_w(filename + "." + self.format)
            title = raw_input("Enter title for latex output: ")
            file.write( header(title,'tex') )
            file.write( self.root.write( \
                  self.format,self.namespace,self.access,recursive=1) )
            file.write( footer('tex') )
            file.close()

        # Copy css stylesheet or xml dtd
        if self.format == 'html' or self.format == 'xml' :
            if not self.src_dir == self.doc_dir:
                ext = ''
                if self.format == 'html' :
                    ext = '.css'
                elif self.format == 'xml' :
                    ext = '.dtd'
                if ext :
                    src  = GDoc_dir     + os.sep + "GDoc" + ext
                    dest = self.doc_dir + os.sep + "GDoc" + ext
                shutil.copy(src,dest)
            
    #end Document.write

#end class Document


# Regular expression pattern for blank line
_blank_line = re.compile(r"^\s*$")


#****c GDoc/InputFile --------------------------------------------------
# CLASS
#     InputFile
# INSTANCE ATTRIBUTES
#     file - file object
#     line - last line read from file, without newline
# METHODS
#     __init__(filename)- open file filename, assign to self.file 
#     readline()  - read and return next line, store in self.line
#     next_line() - find, store, and return the next non-blank line
#     read_string_dict() - read a dictionary, one key value pair per line
#     read_string_list() - read a string list, one string per line
#*** ------------------------------------------------------------------
class InputFile:

    def __init__(self,filename):
        self.file = open(filename,'r')

    #****m InputFile/InputFile.readline  ------------------------------
    # METHOD
    #    InputFile.readline()
    # PURPOSE
    #    Reads and returns one line of file, saves in self.line
    #*** --------------------------------------------------------------
    def readline(self):
        self.line = self.file.readline()
        return self.line

    #****m InputFile/InputFile.next_line  -----------------------------
    # METHOD
    #    InputFile.next_line()
    # PURPOSE
    #    Returns next non-blank line
    #***  -------------------------------------------------------------
    def next_line(self):
        global _blank_line
        self.readline()
        while 1:
            if _blank_line.match(self.line):
                self.readline()
            else:
                return string.strip(self.line)

    #****m InputFile/InputFile.read_string_dict  ----------------------
    # METHOD
    #    InputFile.read_string_dict()
    # PURPOSE
    #    Reads a dictionary of key and value string pairs, from a
    #    file containing one key / value pair per line, with strings
    #    delimited by blank space, until a blank line is encountered.
    #    Returns a dictionary of key and value strings. 
    #***  -------------------------------------------------------------
    def read_string_dict(self):
        string_dict = {}
        self.line = self.next_line()
        while 1:
            if not _blank_line.match(self.line):
                words = string.split(self.line)
                key   = words[0] 
                value = words[1] 
                string_dict[key] = value
                self.line  = self.file.readline()
            else:
                break
        return string_dict


    #****m InputFile/InputFile.read_string_list  -----------------------
    # METHOD
    #    InputFile.read_string_dict()
    # PURPOSE
    #    Reads and a list of strings, with one string per line,
    #    until a blank line is returned. Returns a list of strings
    #    with leading and trailing white space removed.
    #***  -------------------------------------------------------------
    def read_string_list(self) :
        string_list = []
        self.line = self.next_line()
        while 1:
            if not _blank_line.match(self.line):
                self.line = string.strip(self.line)
                string_list.append(self.line)
                self.line = self.file.readline()
            else:
                break
        return string_list

    #****m InputFile/InputFile.read_path -----------------------------
    # METHOD
    #    InputFile.read_path([root])
    # PURPOSE
    #    Read absolute or relative path from next line of file, and 
    #    and return a corresponding absolute path.  If input path name 
    #    is absolute, return the input. Otherwise, treat input path as
    #    a relative path, relative to argument root (if present) or to 
    #    the current working directory (if root is not present).
    #***  -------------------------------------------------------------
    def read_path(self, root='') :
        path = self.next_line()
        if os.path.isabs(path):
            return path
        else:
            if not root: root = os.getcwd()
            return os.path.join(root,path)
            
#end Class InputFile


#****f GDoc/read_config ------------------------------------------------
# FUNCTION
#     read_config(rc_filename)
# PURPOSE
#     Reads data from file rc_filename.
#     Return a Document object with specified values for attributes:
#        src_dir, doc_dir, filenames, format, and access
#     Also sets values for Item and Section class attributes:
#        Item.types, Section.types, Section.ignore
# FORMAT
#    The configuration file consists of a series of control lines,
#    most of which are followed by one or more lines containing data. 
#    There may be no blank lines between the control line and the 
#    associated data line(s), or between data lines. One or more blank
#    lines separate each control line and associated data line(s) from 
#    the next control line.  The recognized control lines and associated 
#    data are:
#
#     'Item_types:'    - dictionary of recognized Item types and 
#                        associated one character codes, e.g.,
#     'Section_types:' - list of recognized section headings
#     'Com_char:'      - character used at beginning of item marker lines
#     'Src_dir:'       - source directory path
#     'Doc_dir:'       - documentation directory path
#     'Filenames:'     - list of names of source files in Src_dir
#     'Format:'        - format of output. Either 'html' or 'xml'
#     'Access:'        - if 'private', include internal Items
#                      - if 'public', exclude internal Items
#     'Done'           - stop reading file, return from function
#
#   Lists (e.g., Section_types, and filenames) or dictionaries
#   (e.g., Item_types) are listed with one item, or one key/value pair,
#   per line.  Directory paths may be either absolute or relative to 
#   current working directory.  The configuration file is read as an 
#   infinite loop, which exits when the "Done" control line is read.
#
#   For examples of the required format, see the GDoc_rc.html or 
#   GDoc_rc.latex file in the GDoc directory, which are used to 
#   generate the html and latex versions of the documentation for
#   the GDoc module itself. 
#*** ------------------------------------------------------------------
def read_config(rc_filename):
    rc_file = InputFile(rc_filename)
    top_dir = os.getcwd()

    # Set defaults
    format = 'html'
    access = 'public'

    # Loop over control lines in rc_file
    while 1:
	line = rc_file.next_line()
        if line == 'Item_types:' :
            Item.types = rc_file.read_string_dict()
            type_patt = ''
            for key in Item.types.keys() :
                type_patt = type_patt + key
            type_patt = r'[' + type_patt + r']'
        elif line == 'Section_types:' :
            Section.types = rc_file.read_string_list()
        elif line == 'Section_ignore:' :
            Section.ignore = rc_file.read_string_list()
        elif line == 'Src_dir:' :
            src_dir = rc_file.read_path()
            print src_dir
        elif line == 'Doc_dir:' :
            doc_dir = rc_file.read_path()
            print doc_dir
        elif line == 'Filenames:' :
            filenames  = rc_file.read_string_list()
        elif line == 'Format:':
            format = rc_file.next_line()
        elif line == 'Access:' :
            access = rc_file.next_line()
        elif line == 'Com_char:' :
            com_char = rc_file.next_line()
        elif line == 'Done':
            # Make comment and Item begin and end marker regex's
            # Item.com, Item.bgn, Item.end = make_markers(com_char,Item.types)
            Item.com = r"^\s*" + com_char
            Item.bgn =  Item.com + r'\*\*\*\*i?' + type_patt
            Item.end =  Item.com + r'\*\*\*'
            Item.com = re.compile(Item.com)
            Item.bgn = re.compile(Item.bgn)
            Item.end = re.compile(Item.end)
            # Construct and return Document object
            Doc = Document(src_dir,doc_dir,filenames,format,access)
            return Doc 
        else :
            raise "Error in read_config: \n " + \
                  "Invalid control line in configuration file:" + \
                   line
#end function read_config

# ------------------- Utility Functions ------------------------------ #

#****f GDoc/link_namespace ---------------------------------------------
# FUNCTION
#    link_namespace(line,namespace)
# ARGUMENTS
#    line      - a string
#    namespace - dictionary of Items accessed by name
# RETURN
#    string line with Item names replaced by href links
#*** -------------------------------------------------------------------
def link_namespace(line,namespace):
    for name in namespace.keys():
        an_Item = namespace[name]
        line = an_Item.html_link(line)
    return line
#end function link_namespace


#****f GDoc/escape ---------------------------------------------------
# FUNCTION
#     escape(a_string,format)
# ARGUMENTS
#     a_string - string
#     format   - output format. May be 'html', 'xml', or 'latex'
# RETURN
#     Returns a copy of a_string appropriate for output in format, 
#     in which special characters are replaced by escape sequences.
#*** -----------------------------------------------------------------
def escape(a_string,format):
    if format == 'html' :
        a_string = a_string.replace(r'&',r'&amp;')
        a_string = a_string.replace(r'<',r'&lt;')
        a_string = a_string.replace(r'>',r'&gt;')
    elif format == 'tex' :
        a_string = a_string.replace('_',r'\_')
    return a_string


#****f GDoc/header ----------------------------------------------
# FUNCTION
#     header(title,format)
# ARGUMENTS
#     title  - title string for html or latex. Note used in xml.
#     format - output format. May be 'html', 'tex', or 'xml'
# RETURN
#     A string containing a standard header for output file.
#*** -----------------------------------------------------------------
def header(title='',format='html'):
    if not title: title = 'GDoc File'
    s=[]
    if format == 'html' :
        style_path = file_util.relative_path(title,'GDoc.css')
        s.append( r'<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 3.2//EN">' \
                  + '\n' )
        s.append( r"<html>" + '\n' )
        s.append( r"<head>" + '\n' )
        s.append( '<link rel="stylesheet" ' + \
                  'href="' + style_path + '" type="text/css">' + '\n')
        s.append( r"<title>" + title + r"</title>" + '\n' )
        s.append( r'</head>' + '\n' )
        s.append( r'<body bgcolor="#FFFFFF">' + '\n' )
        s.append( r"<br />" + '\n' )
    elif format == 'tex' :
        s.append(r'\documentclass{article}'+'\n')
        #s.append(r'\usepackage{makeidx}'+'\n')
        s.append(r'\oddsidemargin  0.00 in'+'\n')
        s.append(r'\evensidemargin 0.00 in'+'\n')
        #s.append(r'\marginparwidth 1.00 in   '+'\n')
        s.append(r'\topmargin      0.00 in   '+'\n')
        s.append(r'\textwidth      6.50 in'+'\n')
        s.append(r'\textheight     9.00 in'+'\n')
        s.append(r'\setlength{\parindent}{0in}'+'\n')
        s.append(r'\setlength{\parskip}{.02in}'+'\n')
        #s.append(r'\pagestyle{headings}'+'\n')
        #s.append(r'\makeindex'+'\n')
        s.append(r'\begin{document}'+'\n')
        s.append(r'\title{' + title + '}' + '\n' )
        s.append(r'\maketitle'+'\n')
        #s.append(r'\printindex'+'\n')
        s.append(r'\tableofcontents'+'\n')
        s.append(r'\newpage'+ '\n')
    elif format == 'xml' :
        s.append( r'<?xml version="1.0" encoding="UTF-8"?>'+'\n')
        s.append( r'<!DOCTYPE SrcFile SYSTEM "GDoc.dtd" >'+'\n')
    return string.join(s,'')


#****f GDoc/footer --------------------------------------------------
# FUNCTION
#     footer(format)
# ARGUMENTS
#     format - output format. 
# PURPOSE
#     Returns a string containing a standard footer for output
#     files in the specified format. Outputs null string for 'xml'
#*** ----------------------------------------------------------------
def footer(format='html'):
    s=[]
    if format == 'html' :
        s.append( r"</body>" + '\n' )
        s.append( r"</html>"  )
    elif format == 'tex' :
        s.append( r'\end{document}' )
    return string.join(s,'')


#****f GDoc/find_directory --------------------------------------
# FUNCTION
#    find_directory(name)
# PURPOSE
#    Searches python path and subdirectories of those in this path 
#    for a directory with basename == name, and returns the absolute 
#    path for this directory. Raises an error if the search fails.
#*** ------------------------------------------------------------
def find_directory(name):
    import sys
    # Search directories in python path
    for dir in sys.path:
	dirname, basename = os.path.split(dir)
	if basename == name:
	    return dir
    # Search subdirectories of those in python path
    for dir in sys.path:
        if os.path.isdir(dir):
            for file in os.listdir(dir):
                path = os.path.join(dir,file)
                if os.path.isdir(path):
                    if file == name:
                        return path
    print "Failed to find directory ", name, " in sys.path="
    print sys.path 
    raise "Error: Search failed in function file_directory"


#****s GDoc/Main ------------------------------------------------------
# SCRIPT
#     Main script
# SOURCE

# Find GDoc source directory
GDoc_dir = find_directory('GDoc')

# Find configuration file
rc_filename = raw_input("Enter name of rc configuration file: ")

Doc = read_config(rc_filename)
Doc.read()
Doc.make_indices()
Doc.make_tree()
Doc.write()

#*** ----------------------------------------------------------------
