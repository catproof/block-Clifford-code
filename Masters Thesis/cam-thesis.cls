%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Class ``cam-thesis''
%%
%% Version: v0.2
%% Authors: Jean Martina, Rok Strnisa, Matej Urbas
%% Date: 30/07/2008
%%
%% Copyright (c) 2008-2012, Rok Strniša, Jean Martina, Matej Urbas
%% License: Simplified BSD License
%% License file: ./License
%% Original License URL: http://www.freebsd.org/copyright/freebsd-license.html
%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% General guidelines on which this class is based:
%%
%%     http://www.cl.cam.ac.uk/local/phd/typography/
%%     http://www.admin.cam.ac.uk/offices/gradstud/exams/submission/phd/format.html
%%
%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Class identification.
%%
%%%%%
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{cam-thesis}[2012/04/12 University of Cambridge thesis class]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Package options (see README.md for a list of options with descriptions).
%%
%% These options can be provided within square brackets of the `documentclass'
%% command.
%%
%%%%%

% times - tells the class to use the times font.
\newif\ifcam@times\cam@timesfalse
\DeclareOption{times}{\cam@timestrue}

% listoffigures - puts the list of figures (after the TOC).
%
\newif\ifcam@listoffigures\cam@listoffiguresfalse
\DeclareOption{listoffigures}{\cam@listoffigurestrue}

% listoftables - puts the list of tables (after the TOC).
%
\newif\ifcam@listoftables\cam@listoftablesfalse
\DeclareOption{listoftables}{\cam@listoftablestrue}

% glossary - puts the glossary (after the TOC).
%
\newif\ifcam@glossary\cam@glossaryfalse
\DeclareOption{glossary}{\cam@glossarytrue}

% index - puts the index at the end of the thesis.
%
\newif\ifcam@index\cam@indexfalse
\DeclareOption{withindex}{\cam@indextrue}

% final - puts the index at the end of the thesis.
%
\newif\ifcam@final\cam@finalfalse
\DeclareOption{final}{\cam@finaltrue}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Using report class as base.
%%
%%%%%
\PassOptionsToClass{12pt}{report}
\ifcam@final
    \PassOptionsToClass{final}{report}
\fi
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\ProcessOptions\relax
\LoadClass{report}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% For consistent vertical spacing
%%
%%%%%
\raggedbottom


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Additional packages, and their options.
%%
%%%%%
\RequirePackage[utf8]{inputenc}
\RequirePackage{calc}
\RequirePackage{microtype}          % Improvement of general text appearance
\RequirePackage[nottoc]{tocbibind}

% Hyperref should be included in any case, but we make it non-obtrusive (i.e.: no hyperlink borders).

\RequirePackage[
    pdffitwindow=true,
    pdfpagelabels=true,
    colorlinks=false,
    pdfborder={0 0 0},
    pdfusetitle=true
]{hyperref}
\RequirePackage[all]{hypcap}        % Fixes links to floats (must be included after hyperref).

% Create the glossary
\ifcam@glossary
    \RequirePackage[toc,nonumberlist,acronyms]{glossaries}
    \makeglossaries%
    \setglossarystyle{listdotted}
\fi


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Page margins.
%%
%%%%%

\RequirePackage{geometry}
\newlength{\cam@bottom}
\newlength{\cam@marginparwidth}
\let\oldgeometry\geometry
\let\oldnewgeometry\newgeometry
\renewcommand{\geometry}[5][0.7]{
    \setlength{\cam@marginparwidth}{#2}
    \addtolength{\cam@marginparwidth}{-2.5mm}
    \setlength{\cam@bottom}{#5}
    \oldgeometry{letterpaper,left=#2,right=#3,top=#4,bottom=\cam@bottom+#1\cam@bottom,footskip=#1\cam@bottom,marginparwidth=\cam@marginparwidth,marginparsep=2mm}
}
\renewcommand{\newgeometry}[5][0.7]{
    \setlength{\cam@marginparwidth}{#2}
    \addtolength{\cam@marginparwidth}{-2.5mm}
    \setlength{\cam@bottom}{#5}
    \oldnewgeometry{left=#2,right=#3,top=#4,bottom=\cam@bottom+#1\cam@bottom,footskip=#1\cam@bottom,marginparwidth=\cam@marginparwidth,marginparsep=2mm}
}
\geometry{1.5in}{1in}{1in}{1in}
\reversemarginpar


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Line spacing
%%
%%%%%

\RequirePackage[doublespacing]{setspace}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Todos
%%
%%%%%

\RequirePackage[obeyFinal,colorinlistoftodos]{todonotes}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Variable definitions and default values: these variables should be defined by
%% the user (somewhere in the preamble). For example, to put the abstract into
%% the thesis, the thesis writer should type the following somewhere in the
%% preamble (before the `\begin{document}` or `\frontmatter` commands are
%% called):
%%
%%      \abstract{This is my abstract.}
%%
%% See below (in the comments starting with 'DOCVAR: ') for a list of all
% variables
%% the thesis writer is expected to use.
%%
%%%%%

% DOCVAR: abstract (The text that will be instered into the abstract of the
% thesis.)
\newcommand{\@abstract}{}
\renewcommand{\abstract}[1]{\renewcommand{\@abstract}{#1}}

% DOCVAR: acknowledgements (The text that will be instered into the
% acknowledgments of the thesis.)
\newcommand{\@acknowledgements}{}
\newcommand{\acknowledgements}[1]{\renewcommand{\@acknowledgements}{#1}}

% DOCVAR: college (The name of the thesis writer's college, which will appear
% just below their name.)
\newcommand{\@college}{}
\newcommand{\college}[1]{\renewcommand{\@college}{#1}}

% DOCVAR: location (The location of the thesis writer's college, which will appear
% just below their name.)
\newcommand{\@location}{}
\newcommand{\location}[1]{\renewcommand{\@location}{#1}}

% DOCVAR: keywords (These keywords will appear in the PDF meta-information
% called `pdfkeywords`.)
\newcommand{\@keywords}{}
\newcommand{\keywords}[1]{\renewcommand{\@keywords}{#1}}

% DOCVAR: subjectline (This subject will appear in the PDF meta-information
% called `pdfsubject`.)
\newcommand{\@subjectline}{}
\newcommand{\subjectline}[1]{\renewcommand{\@subjectline}{#1}}

% DOCVAR: submissiondate (The date of the submission of this thesis.)
\newif\ifcam@submissiondate\cam@submissiondatefalse
\newcommand{\@submissiondate}{}
\newcommand{\submissiondate}[1]{\renewcommand{\@submissiondate}{#1}\cam@submissiondatetrue}

% DOCVAR: type (The type of document, e.g., thesis, thesis proposal, dissertation.)
\newcommand{\@thesistype}{}
\newcommand{\thesistype}[1]{\renewcommand{\@thesistype}{#1}}

% DOCVAR: type (The type of document, e.g., thesis, thesis proposal, dissertation.)
\newcommand{\@submittedto}{}
\newcommand{\submittedto}[1]{\renewcommand{\@submittedto}{#1}}

% DOCVAR: submissionnotice (The submission notice is shown on the bottom of the
% title page.)
% Faculty of Graduate and Postdoctoral Affairs
\newcommand{\@submissionnotice}{%
A {\@thesistype} submitted to {\@submittedto} in partial fulfillment of the requirements for the degree of
}
\newcommand{\submissionnotice}[1]{\renewcommand{\@submissionnotice}{#1}}

% DOCVAR: degree (The degree for which this thesis is written.)
\newcommand{\@degree}{}
\newcommand{\degree}[1]{\renewcommand{\@degree}{#1}}

% DOCVAR: program (The program for which this thesis is written.)
\newcommand{\@program}{}
\newcommand{\program}[1]{\renewcommand{\@program}{#1}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Chapter and section numbering
%%
% TODO: Check whether numbering up to the fourth level is okay. It is not
% specified in the guidelines.
\setcounter{secnumdepth}{2}
\setcounter{tocdepth}{2}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Front matter
%%
%% - outside and inside front cover
%% - title leaf
%% Do not include the date of make!
%% Institution + department.
%% Names of referees. (optional)
%% Degree.
%% Date of submission and defense. (optional)
%% Place and date of publication and publishers (and other info by them).
%%%%%
\newcommand{\frontmatter}{

\pagestyle{plain}
\newgeometry[0]{1.5in}{1.5in}{1.5in}{1.5in}

\ifcam@final
\else
    \pagenumbering{roman}
    \setcounter{page}{0}
    \thispagestyle{empty}
    \listoftodos
    \newpage
\fi

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Title Page
%%
%%%%%
\pagenumbering{roman}
\setcounter{page}{0}
\thispagestyle{empty}
\begin{center}
    {
        \Huge\bfseries
        \@title
    }
    \vfill
    {
        by
    }
    \vfill
    {
        \Large\bfseries
        \@author
    }
    \vfill
    {
        \@submissionnotice
    }
    \vfill
    {
        \large\bfseries
        \@degree
    }
    \vfill
    {
        in
    }
    \vfill
    {
        \large\bfseries
        \@program
    }
    \vfill
    {
        \@college\\
        \@location
    }
    \vfill
    {
        \textcopyright~\@submissiondate\\
        \@author
    }
\end{center}
\hypersetup{pdfsubject={\@subjectline},pdfkeywords={\@keywords}}

\newpage
\restoregeometry

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Abstract
%%
%%%%%

\chapter*{Abstract}
\addcontentsline{toc}{chapter}{Abstract}
\@abstract{}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Acknowledgements
%%
%%%%%

\chapter*{Acknowledgements}
\addcontentsline{toc}{chapter}{Acknowledgements}
\@acknowledgements{}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Table of contents, figures, symbols and glossary.
%%
%%%%%

\tableofcontents

\ifcam@listoftables%
    \listoftables
\fi

\ifcam@listoffigures%
    \listoffigures
\fi

\ifcam@glossary%
    \printglossaries
\fi

\newpage
\setcounter{page}{1}
\pagenumbering{arabic}
}
