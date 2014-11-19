NAME = crude1
CC = gcc
CFLAGS = -g
MLIB = -lm

$(NAME): $(NAME).c
	$(CC) $(C_FLAGS) -o $(NAME) $(NAME).c $(MLIB)
c: $(NAME).w
	ctangle $(NAME).w
tex: $(NAME).w
	cweave $(NAME).w
dvi: $(NAME).tex
	tex $(NAME).tex
ps: $(NAME).dvi
	dvips $(NAME).dvi
pdf: $(NAME).dvi
	dvipdf $(NAME).dvi
clean:
	rm -f $(NAME).c $(NAME).tex $(NAME).dvi $(NAME).ps $(NAME).pdf $(NAME).idx $(NAME).scn $(NAME).toc $(NAME).log $(NAME).*~ $(NAME)


