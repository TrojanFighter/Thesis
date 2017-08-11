%%
%% This is file `cufethesis.cls',
%% generated with the docstrip utility.
%%
%% The original source files were:
%%
%% cufethesis.dtx  (with options: `class')
%%
%% Copyright(C) 2014 by Mohsen Muhammad Mahroos, PhD
%%
%% Assistant Professor
%% Electronics and Electrical Communications Engineering Department
%% Faculty of Engineering
%% Cairo University
%% Giza, EGYPT 12613
%% Phone: +202 3567-8888
%% E-mail: mmahroos@eece.cu.edu.eg
%%
%% This is a free file may be distributed and/or modified under the
%% conditions of the latest version of the LaTeX Project Public License
%% available in:
%%
%%           http://www.latex-project.org/lppl.txt
%%
%% which is part of all distributions of LaTeX on June 17, 2014.
%%
%% This is a generated file designed to implement LaTeX Document Class for
%% typesetting MSc and PhD Theses prepared in all Departments in the Faculty
%% of Engineering at Cairo University according to the Graduate Studies
%% Administration Guidelines released in January 2013.
%%
\NeedsTeXFormat{LaTeX2e}[2011/06/27]
\ProvidesClass{cufethesis}
  [2014/06/17 v1.2 A LaTeX Document Class for Typesetting Graduate Theses ]
%% \CharacterTable
%%  {Upper-case    \A\B\C\D\E\F\G\H\I\J\K\L\M\N\O\P\Q\R\S\T\U\V\W\X\Y\Z
%%   Lower-case    \a\b\c\d\e\f\g\h\i\j\k\l\m\n\o\p\q\r\s\t\u\v\w\x\y\z
%%   Digits        \0\1\2\3\4\5\6\7\8\9
%%   Exclamation   \!     Double quote  \"     Hash (number) \#
%%   Dollar        \$     Percent       \%     Ampersand     \&
%%   Acute accent  \'     Left paren    \(     Right paren   \)
%%   Asterisk      \*     Plus          \+     Comma         \,
%%   Minus         \-     Point         \.     Solidus       \/
%%   Colon         \:     Semicolon     \;     Less than     \<
%%   Equals        \=     Greater than  \>     Question mark \?
%%   Commercial at \@     Left bracket  \[     Backslash     \\
%%   Right bracket \]     Circumflex    \^     Underscore    \_
%%   Grave accent  \`     Left brace    \{     Vertical bar  \|
%%   Right brace   \}     Tilde         \~}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{report}}
\def\@MSc{Master of Science}
\DeclareOption{msc}{\def\@Degree{\@MSc}}
\def\@PhD{Doctor of Philosophy}
\DeclareOption{phd}{\def\@Degree{\@PhD}}
\def\@Degree{\@MSc}
\DeclareOption{emp}{
    \def\@Specialization{Engineering Mathematics and Physics}
}
\DeclareOption{arch}{
    \def\@Specialization{Architectural Engineering}
}
\DeclareOption{civil}{
    \def\@Specialization{Civil Engineering}
}
\DeclareOption{mech}{
    \def\@Specialization{Mechanical Engineering}
}
\DeclareOption{aero}{
    \def\@Specialization{Aerospace Engineering}                               % modified by mhawwary
}
\DeclareOption{eece}{
    \def\@Specialization{Electronics and Communications}
}
\DeclareOption{epm}{
    \def\@Specialization{Electrical Power and Machines}
}
\DeclareOption{chem}{
    \def\@Specialization{Chemical Engineering}
}
\DeclareOption{mpm}{
    \def\@Specialization{Mining, Petroleum and Metallurgical}
}
\DeclareOption{bio}{
    \def\@Specialization{Biomedical Engineering and Systems}
}
\DeclareOption{comp}{
    \def\@Specialization{Computer Engineering}
}
\def\@Specialization{Aerospace Engineering}

\newif\if@praise
\DeclareOption{nopraise}{\@praisefalse}
\DeclareOption{praise}{\@praisetrue}
\newif\if@phero
\DeclareOption{nophero}{\@pherofalse}
\DeclareOption{phero}{\@pherotrue}
\newif\if@longtitle
\DeclareOption{nolongtitle}{\@longtitlefalse}
\DeclareOption{longtitle}{\@longtitletrue}
\newif\if@infopage
\DeclareOption{noinfopage}{\@infopagefalse}
\DeclareOption{infopage}{\@infopagetrue}
\newif\if@listoftables
\DeclareOption{nolistoftables}{\@listoftablesfalse}
\DeclareOption{listoftables}{\@listoftablestrue}
\newif\if@listoffigures
\DeclareOption{nolistoffigures}{\@listoffiguresfalse}
\DeclareOption{listoffigures}{\@listoffigurestrue}
\newif\if@listofalgorithms
\DeclareOption{nolistofalgorithms}{\@listofalgorithmsfalse}
\DeclareOption{listofalgorithms}{\@listofalgorithmstrue}

%----------------------------------------% added by mhawwary-------------------------------------------------------
\newif\if@acknow
\DeclareOption{noacknow}{\@acknowfalse}
\DeclareOption{acknow}{\@acknowtrue}

\newif\if@dedic
\DeclareOption{nodedic}{\@dedicfalse}
\DeclareOption{dedic}{\@dedictrue}

%------------------------------------------------------------------------------------------------------------------------------

% \ExecuteOptions{praise,nophero,nolongtitle,infopage}     
\ExecuteOptions{praise,acknow,dedic,nophero,nolongtitle,infopage}    % modified by mhawwary
\ExecuteOptions{listoftables,listoffigures,nolistofalgorithms}
\ProcessOptions

\LoadClass[12pt,a4paper,twoside]{report}
\setlength\oddsidemargin{10mm}
\setlength\evensidemargin{0mm}
\setlength\topmargin{0mm}
\setlength\headheight{0mm}
\setlength\textheight{238.2mm}
\setlength\textwidth{149.2mm}
\setcounter{tocdepth}{4}
\setcounter{secnumdepth}{4}
\lefthyphenmin 62
\righthyphenmin 62

%\usepackage[inner=10mm,outer=5mm]{geometry}
\usepackage[top = 1.0in, bottom = 1.25in, left = 1.3in, right = 0.9in]{geometry}      % modified by mhawwary

\usepackage{alltt}
\usepackage[titlenotnumbered,boxruled]{algorithm2e}
\usepackage[centertags]{amsmath}
\usepackage{amsfonts}
\usepackage{amssymb}
\usepackage{amsthm}
\usepackage[toc,page]{appendix}
%\usepackage[titletoc]{appendix}
\usepackage[english]{babel}
\usepackage{calc}
\usepackage{epsfig}
\usepackage{psfrag}
\usepackage[breaklinks]{hyperref}
%\hypersetup{   linktoc=all  }
\usepackage[all]{hypcap}
\usepackage{ifthen}
\usepackage{pdflscape}
\usepackage{mathptmx}
\usepackage{microtype}
\usepackage{multicol}
\usepackage{multirow}
\usepackage{newlfont}
\usepackage[intoc]{nomencl}
\usepackage[overload]{textcase}
%\usepackage{tocloft}
\usepackage[subfigure]{tocloft}
\usepackage{txfonts}
\usepackage[notbib]{tocbibind}
\usepackage{xstring}
\usepackage{stringstrings}
\usepackage{graphicx}
\usepackage{epstopdf}
\usepackage{wrapfig}
 \usepackage{subfigure}% subcaptions for subfigures
 %\usepackage[subfigure]{tocloft}
 
 \usepackage{subfigmat}% matrices of similar subfigures, aka small mulitples
 \usepackage{varioref}%  smart page, figure, table, and equation referencing
\usepackage{float}
\usepackage{algorithm2e}
%\usepackage{hyperref}
%\usepackage{caption}
\usepackage{latexsym}
%\usepackage{mathrsfs}
%\usepackage{natbib}
%\usepackage{textcomp}
\usepackage{pifont}
\usepackage{varioref}
\usepackage{color}
\usepackage[font=normalsize,skip=2.0pt]{caption}
\usepackage{arabtex}
\usepackage{setspace}     \singlespacing
\usepackage{tocloft}
%\usepackage{etoolbox}

\graphicspath{{fig/}}
\graphicspath{{../fig/}}


 \newcommand{\BibTeX}{\textsc{Bib}\TeX}

\def\specialization{\@Specialization}
\def\department{\@Specialization~Department}
\def\faculty{Faculty of Engineering}
\def\university{Cairo University}
\def\city{Giza}
\def\country{Egypt}
\newcommand\MainSupervisor[2]
{
    \def\@MainSupervisor{#1}
    \def\@MainSupervisorPosition{#2}
}
\newcounter{@NumberOfCoSupervisors}

\newcommand{\CoSupervisor}[4]{
    \addtocounter{@NumberOfCoSupervisors}{1}
    \ifcase\the@NumberOfCoSupervisors\relax
    % Do nothing if no Co-Supervisors
    \or \def\@FirstCoSupervisor{#1}
        \def\@FirstCoSupervisorPosition{#2}
        \def\@FirstCoSupervisorDepartment{#3}
        \def\@FirstCoSupervisorAffiliation{#4}
    \or \def\@SecondCoSupervisor{#1}
        \def\@SecondCoSupervisorPosition{#2}
        \def\@SecondCoSupervisorDepartment{#3}
        \def\@SecondCoSupervisorAffiliation{#4}
    \or \def\@ThirdCoSupervisor{#1}
        \def\@ThirdCoSupervisorPosition{#2}
        \def\@ThirdCoSupervisorDepartment{#3}
        \def\@ThirdCoSupervisorAffiliation{#4}
    \fi
}
\newcommand{\ExternalExaminer}[1]{\def\@ExternalExaminer{#1}}
\newcommand\InternalExaminer[1]{\def\@InternalExaminer{#1}}
\newcounter{@NumberOfMemberExaminers}

\newcommand\MemberExaminer[1]{
    \addtocounter{@NumberOfMemberExaminers}{1}
    \ifcase\the@NumberOfMemberExaminers\relax
    % Do nothing if no Member Examiners
    \or \def\@FirstMemberExaminer{#1}
    \or \def\@SecondMemberExaminer{#1}
    \fi
}
\newcommand{\RegistrationDate}[3]{\def\@RegistrationDate{#1/#2/#3}}
\newcommand{\SubmissionDate}[2]{
    \def\@SMonth{#1}
    \def\@SYear{#2}
     %\def\@SYear{#1}
}
\newcommand{\AwardingDate}[3]{\def\@AwardingDate{#1/#2/#3}}
\def\@AwardingDate{}
\newcommand{\EngineerName}[1]{\def\@EngineerName{#1}}
\newcommand{\BirthDate}[3]{\def\@BirthDate{#1/#2/#3}}
\newcommand{\Nationality}[1]{\def\@Nationality{#1}}
\newcommand{\EMail}[1]{\def\@EMail{#1}}
\newcommand{\Phone}[1]{\def\@Phone{#1}}
\newcommand{\Address}[1]{\def\@Address{#1}}
\newcommand{\Photo}[1]{\def\@Photo{#1}}
\newcommand{\ThesisTitle}[1]{\def\@ThesisTitle{#1}}
\newcommand{\ThesisKeywords}[1]{\def\@ThesisKeywords{#1}}
\newcommand{\ThesisSummary}[1]{\def\@ThesisSummary{#1}}


\setlength{\cftbeforesubsubsecskip}{0.05in}
\setlength{\cftbeforesubsecskip}{0.1in}
\setlength{\cftbeforesecskip}{0.15in}
\setlength{\cftbeforechapskip}{0.275in}

\usepackage[svgnames]{xcolor}
\makeatletter
\let\stdl@section\l@section
\renewcommand*{\l@section}[2]{%
  \stdl@section{\uppercase{#1}}{\textcolor{black}{#2}}}

\newcommand{\Preamble}[1]{
    \@NNChapter{#1}
    \input{#1.tex}
    \newpage
}
\newif\if@firstchapter\@firstchaptertrue


\newcommand{\Chapter}[2]{ %
    \gchapter{\caseupper{#1}}
    \chapter{#1}
    \label{Chapter:#2}
    \if@firstchapter
        \@firstchapterfalse
        \pagestyle{plain}
        \pagenumbering{arabic}
    \fi
    \input{#2.tex}
    \newpage
}


\newcommand{\Appendix}[2]{
    \gappendix{\caseupper{#1}}
    \chapter{#1}
    \label{Appendix:#2}
    \input{#2.tex}
    \newpage
}

\renewcommand\chapter{\if@openright\cleardoublepage\else\clearpage\fi
                    \thispagestyle{plain}%
                    \global\@topnum\z@
                    \@afterindentfalse
                    \secdef\@chapter\@schapter }
\def\@chapter[#1]#2{\chaptermark{#1}%
                    \addtocontents{lof}{\protect\addvspace{10\p@}}%
                    \addtocontents{lot}{\protect\addvspace{10\p@}}%
                    \if@twocolumn
                      \@topnewpage[\@makechapterhead{#2}]% 
                    \else
                      \@makechapterhead{#2}%
                      \@afterheading
                    \fi}

%------------------------------------------------------------------------------------
%                            Labeling and Referencing Commands
%-----------------------------------------------------------------------------------

\newcommand{\seclbl}[1]{\label{Section:#1}}

\newcommand{\chref}[1]{Chapter~\ref{Chapter:#1}}
\newcommand{\appref}[1]{Appendix~\ref{Appendix:#1}}

\newcommand{\secref}[1]{Section~\ref{Sec:#1}}
\newcommand{\eqnlbl}[1]{\label{Equation:#1}}
\newcommand{\eqnref}[1]{Equation~\ref{Eqn:#1}}
\newcommand{\tabref}[1]{Table~\ref{Table:#1}}
\newcommand{\figref}[1]{Figure~\ref{Fig:#1}}
\newcommand{\algref}[1]{Algorithm~\ref{Algorithm:#1}}

%--------------------------------------------------------------------------------

\newcommand{\gtable}[4]{
    \begin{hctable}{#1}{#2}
        \includegraphics[#3]{#4}
    \end{hctable}
}



\newcommand{\gfigure}[4]
{
    \begin{hcfigure}
\centering
        \includegraphics[#3]{#4}
\centering
        \caption{#1}
        \label{Figure:#2}
    \end{hcfigure}
}
\newenvironment{Algorithm}[2]{
    \begin{algorithm}[H]
        \caption{#1}
        \label{Algorithm:#2}
}
{
    \end{algorithm}
}
\newcommand{\newnom}[2]{
    \nomenclature{#2}{#1}#1~(#2)
}


%------------------------------------------------------------------------------------
%                            Thesis Title Page Command
%-----------------------------------------------------------------------------------
\newenvironment{@thesistitlepage}{
    \thispagestyle{empty}%
    \begin{center}
    \if@longtitle
      \large
    \else
      \Large
    \fi
    \textbf{\MakeTextUppercase{\@ThesisTitle}}\\[28pt]
    \Large
    \textrm{By}\\[16pt]
    \textrm{\@EngineerName}\\[28pt]
    \large A Thesis Submitted to the\\
        \faculty \space at \university\\
        in Partial Fulfillment of the\\
        Requirements for the Degree of\\
        \MakeTextUppercase{\@Degree}\\
        in \\ \@Specialization \\[28pt]
}
{   \vfill\large
    \centerline{\MakeTextUppercase{\department}}
    \centerline{\MakeTextUppercase{\faculty},
                \MakeTextUppercase{\university}}
    \centerline{\MakeTextUppercase{\city},
                \MakeTextUppercase{\country}} 
    %\centerline{\MakeTextUppercase{\@SMonthText} {\@SYear}}
     \centerline{\MakeTextUppercase {\@SYear}}
    \normalsize
    \end{center}
    \newpage
}

%------------------------------------------------------------------------------------
%                            Cover Page Command
%-----------------------------------------------------------------------------------
\newcommand{\@MakeEnglishCoverPage}{
    \if@praise
        \centerline{\includegraphics{Praise}}
    \fi
    \includegraphics{fig/logo/FE-Logo} \hfill
    \if@phero
        \includegraphics{fig/logo/CU-Phero}
    \else
        \includegraphics{fig/logo/CU-Logo}
    \fi
    \vspace*{1.0in}
    \begin{@thesistitlepage}
    \end{@thesistitlepage}
}

%------------------------------------------------------------------------------------
%                            Supervision Committee Page Command
%-----------------------------------------------------------------------------------
\newcommand{\@MakeEnglishSupervisionCommitteePage}{
    \begin{@thesistitlepage}
        Under the Supervision of\\[14pt] 
        \ifcase\the@NumberOfCoSupervisors
            % No Co-Supervisors
            \begin{tabular}{c}
                \@MainSupervisor\\ ....................... \\
                \small\@MainSupervisorPosition\\
                \small\department\\
                \small\faculty, \university\\
            \end{tabular}
        \or % One Co-Supervisor
            \begin{tabular}{cc}
                \@MainSupervisor &
                \@FirstCoSupervisor\\ ........................ &  ....................... \\
                \footnotesize\@MainSupervisorPosition &
                \footnotesize\@FirstCoSupervisorPosition\\
                \footnotesize\department&
                \footnotesize\@FirstCoSupervisorDepartment \\
                \footnotesize\faculty, \university &
                \footnotesize\@FirstCoSupervisorAffiliation\\
            \end{tabular}
        \or % Two Co-Supervisors
            \begin{tabular}{cc}
                \@MainSupervisor &
                \@FirstCoSupervisor\\ ........................ &  ....................... \\
                \footnotesize\@MainSupervisorPosition &
                \footnotesize\@FirstCoSupervisorPosition\\
                \footnotesize\department&
                \footnotesize\@FirstCoSupervisorDepartment\\
                \footnotesize\faculty, \university &
                \footnotesize\@FirstCoSupervisorAffiliation\\
            \end{tabular}
            \begin{tabular}{c} \\
                \@SecondCoSupervisor\\  ......................... \\
                \footnotesize\@SecondCoSupervisorPosition\\
                \footnotesize\@SecondCoSupervisorDepartment\\
                \footnotesize\@SecondCoSupervisorAffiliation
            \end{tabular}
        \or % Three Co-Supervisors
            \begin{tabular}{cc}
                \@MainSupervisor &
                \@FirstCoSupervisor\\
                \footnotesize\@MainSupervisorPosition &
                \footnotesize\@FirstCoSupervisorPosition\\
                \footnotesize\department&
                \footnotesize\@FirstCoSupervisorDepartment\\
                \footnotesize\faculty, \university &
                \footnotesize\@FirstCoSupervisorAffiliation\\
            \end{tabular}
            \begin{tabular}{cc} \\
                \@SecondCoSupervisor &
                \@ThirdCoSupervisor\\
                \footnotesize\@SecondCoSupervisorPosition&
                \footnotesize\@ThirdCoSupervisorPosition\\
                \footnotesize\@SecondCoSupervisorDepartment&
                \footnotesize\@ThirdCoSupervisorDepartment\\
                \footnotesize\@SecondCoSupervisorAffiliation&
                \footnotesize\@ThirdCoSupervisorAffiliation
            \end{tabular}
        \fi
    \end{@thesistitlepage}
}

%------------------------------------------------------------------------------------
%                           Examining Committee Page Command
%-----------------------------------------------------------------------------------
\newcommand{\@MakeEnglishExaminingCommitteePage}{
    \begin{@thesistitlepage}
        \begin{quote}
            Approved by the Examining Committee:\\[24pt]
            \begin{tabular}{l}
                \\  \hline
                \@MainSupervisor, Thesis Main Advisor
                \\ \\ \\
                \ifnum \the@NumberOfMemberExaminers > 0
                    \hline \@FirstMemberExaminer, Member 
                    \\ \\ \\ \hline 
                    \ifnum \the@NumberOfMemberExaminers > 1
                        \hline
                        \@SecondMemberExaminer, Member
                        \\ \\ \\  \hline 
                    \fi
                \fi
 	   \@InternalExaminer, Internal Examiner  \\ \\ \\  \hline
                \@ExternalExaminer, External Examiner
                \\ \\ 
            \end{tabular}
        \end{quote}
    \end{@thesistitlepage}
}

%------------------------------------------------------------------------------------
%                           Acknowledgements Page Command
%-----------------------------------------------------------------------------------
\newcommand{\@MakeAcknowledgementsPage}
{
    \clearpage\phantomsection
    \Preamble{Acknowledgements}
}

%------------------------------------------------------------------------------------
%                           Dedication Page Command
%-----------------------------------------------------------------------------------
\newcommand{\@MakeDedicationPage}
{
    \clearpage\phantomsection
    \Preamble{Dedication}
}

%------------------------------------------------------------------------------------
%                           Abstract Page Command
%-----------------------------------------------------------------------------------
\newcommand{\@MakeAbstractPage}
{
  \clearpage\phantomsection
   \Preamble{Abstract}
}

%------------------------------------------------------------------------------------
%                          English Information Page Command
%-----------------------------------------------------------------------------------
\newcommand{\@MakeEnglishInformationPage}{
    \thispagestyle{empty}
    \begin{quote}
    \begin{tabular}{lll}
        % Personal Information
        \textbf{Engineer's Name:}  & \@EngineerName&
        \multirow{9}{*}{\includegraphics[width=3.4 cm]{\@Photo}}\\
        \textbf{Date of Birth:}     & \@BirthDate&\\
        \textbf{Nationality:}       & \@Nationality&\\
        \textbf{E-mail:}            & \@EMail&\\
        \textbf{Phone:}             & \@Phone&\\
        \textbf{Address:}           & \@Address&\\
        \textbf{Registration Date:} & \@RegistrationDate&\\
        \textbf{Awarding Date:}     & \@AwardingDate&\\
        \textbf{Degree:}            & \@Degree&\\
        \textbf{Department:}        & \@Specialization&\\ & &\\
        \textbf{Supervisors:}       & &\\
                                    & \@MainSupervisor&\\
        \ifnum  \the@NumberOfCoSupervisors > 0
            & \@FirstCoSupervisor&\\
            \ifnum \the@NumberOfCoSupervisors > 1
                & \@SecondCoSupervisor&\\
                    \ifnum \the@NumberOfCoSupervisors > 2
                        & \@ThirdCoSupervisor&\\
                    \fi
            \fi
        \fi
        & &\\
        \textbf{Examiners:}
            & &\\
            & \@ExternalExaminer  &(External examiner)  \\
            & \@InternalExaminer &(Internal examiner) \\
            & \@MainSupervisor  &(Thesis main advisor)\\
            \ifnum \the@NumberOfMemberExaminers > 0
                & \@FirstMemberExaminer &(Member)\\
                \ifnum \the@NumberOfMemberExaminers > 1
                    & \@SecondMemberExaminer   &(Member)\\
                \fi
            \fi
        & & \\
    \end{tabular}

    \textbf{Title of Thesis:}
        \begin{center}
            \@ThesisTitle
        \end{center}

    \textbf{Key Words:}\\ \\ \@ThesisKeywords

    \textbf{Summary:}\\ \@ThesisSummary
    \end{quote}
    \newpage
}

\newcommand{\@SMonthText}{
\ifcase\@SMonth\relax
\or January
\or February
\or March
\or April
\or May
\or June
\or July
\or August
\or September
\or October
\or November
\or December
\fi
}

%-------------------------------------------------------------------
%    Not numbered chapters title's command
%------------------------------------------------------------------------
\newcommand{\@NNChapter}[1]{
    \newpage
    \@NNChapterTitle{#1}
    \addcontentsline{toc}{chapter}{\caseupper{#1}}
}

\newcommand{\@NNChapterTitle}[1]{
    \centerline{\textbf{\Huge #1}}
     \vspace*{24pt}
}

\newcommand{\@NNNChapter}[1]{
    \newpage
    \phantomsection
    \@NNChapterTitle{#1}
    \addcontentsline{toc}{chapter}{\caseupper{#1}}
}
%------------------------------------------------------------------------------

\def\@makechapterhead#1{
  {
  \huge\bfseries
  \settowidth{\@chapterheadwidth}{\@chapapp~\thechapter:~#1}
  \ifnum \@chapterheadwidth < \textwidth
    \begin{center} % The required centering of one-line chapter headers
        \@chapapp~\thechapter:~#1
    \end{center}
  \else % Flush-left multi-line chapter header
    \begin{flushleft}
        \@chapapp~\thechapter: #1
    \end{flushleft}
  \fi
  }
  \vspace*{4pt}
}
\newlength{\@chapterheadwidth}


\newenvironment{hcfigure}
{
     \begin{hfigure}
        \begin{center}
}
{
        \end{center}
     \end{hfigure}
}
\newenvironment{hfigure}{\def\@captype{figure}}{}


\newenvironment{hctable}[2]
{
    \begin{htable}
        \caption{#1}
        \label{Table:#2}
        \begin{center}
            \begin{tabular}{c}
}
{           \end{tabular}
        \end{center}
    \end{htable}
}



\newenvironment{htable}{\def\@captype{table}}{}
\newlength{\nomenclaturelabelwidth}
\setlength{\nomenclaturelabelwidth}{1.0in}
\makenomenclature\makeindex
\def\thenomenclature{%
    \@NNChapter{Nomenclature}
    \nompreamble
    \list{}{%
        \labelwidth\nom@tempdim
        \leftmargin\labelwidth
        \advance\leftmargin\labelsep
        \itemsep\nomitemsep
        \let\makelabel\nomlabel
    }
}


\newcommand{\gchapter}[1]
{
  \ifnum \c@secnumdepth >\m@ne
    \refstepcounter{chapter}
    \typeout{\@chapapp\space\thechapter.}
    \addcontentsline{toc}{chapter}{CHAPTER {\thechapter}: #1}
  \else
    \addcontentsline{toc}{chapter}{#1}
  \fi
}


\newcommand{\gappendix}[1]
{
  \ifnum \c@secnumdepth >\m@ne
    \refstepcounter{chapter}
    \typeout{\@chapapp\space\thechapter.}
    \addcontentsline{toc}{chapter}{APPENDIX {\thechapter}: #1}
  \else
    \addcontentsline{toc}{chapter}{#1}
  \fi
}

\appendixtitletocon\appendixpageoff\appendixtocoff
\numberwithin{equation}{chapter}
\renewcommand{\theequation}{\thechapter.\arabic{equation}}

%------------------------------------------------------------------------------------------------------------------------
%                                                    Begin Document Section 
% -----------------------------------------------------------------------------------------------------------------------
\AtBeginDocument{
    \@MakeEnglishCoverPage                                          % English Cover Page declaration
    \@MakeEnglishSupervisionCommitteePage                % English Supervision Committee Page Declaration
    \@MakeEnglishExaminingCommitteePage                  % English Examining Committee Page Declaration
    \if@infopage
        \@MakeEnglishInformationPage                             % English Examining Committee Page Declaration
    \fi

    \pagenumbering{roman}
    \pagestyle{plain}
    \setcounter{page}{1}
    
    \if@acknow
    \@MakeAcknowledgementsPage              %    Acknowledgements Page declaration, added by mhawwary
    \fi
    \if@dedic
    \@MakeDedicationPage                            %    Dedication Page declaration, added by mhawwary
    \fi
    
    \phantomsection
    \@NNChapter{Table of Contents}             % Table of Contents declaration
    \@starttoc{toc}

    \if@listoftables
        \@NNNChapter{List of Tables}               % List of Tables declaration
        \@starttoc{lot}
    \fi
    \if@listoffigures
        \@NNNChapter{List of Figures}             % List of Figures declaration
        \@starttoc{lof}
    \fi
    \if@listofalgorithms
        \@NNNChapter{List of Algorithms}         % List of Algorithms declaration
        \@starttoc{loa}
    \fi

    \clearpage\phantomsection
    \printnomenclature[\nomenclaturelabelwidth]    % Nomenclature declaration, added by mhawwary

    \@MakeAbstractPage                                          % Abstract Page declaration, added by mhawwary
}
\endinput
%%
%% End of file `cufethesis.cls'.
