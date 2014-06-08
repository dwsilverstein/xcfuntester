      subroutine colorize(text, lerror)      

      implicit none
c
      character*68 text
      character*90 newtext
c
      logical lerror
c
c This routine colorizes text by making use of ANSI escape sequences.
c These allow the terminal to color and bold text.  YMMV.  See list
c below for options other than the bold red and bold black provided
c in this routine (other options probably exist also):
c
c ANSI Code                  Color
c ---------         ----------------------
c    30               Black
c   1;30              Bold Black
c    31               Red
c   1;31              Bold Red
c    32               Green
c   1;32              Bold Green
c    33               Orange/Yellow
c   1;33              Bold Orange/Yellow
c    34               Blue
c   1;34              Bold Blue
c    35               Purple
c   1;35              Bold Purple
c    36               Cyan
c   1;36              Bold Cyan
c
c
      if (lerror) then
        text = adjustl(text)
        text = trim(text)
        write(newtext,*) achar(27) // '[1;31m' // 
     &                   text // achar(27) // '[0m'
        newtext = adjustl(newtext)
        write(*,*) trim(newtext)
      else
        text = adjustl(text)
        text = trim(text)
        write(newtext,*) achar(27) // '[1;30m' //
     &                   text // achar(27) // '[0m'
        newtext = adjustl(newtext)
        write(*,*) trim(newtext)
      endif

      end
