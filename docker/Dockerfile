FROM ubuntu:latest

COPY CADBURE /usr/bin/

RUN apt update && apt install cpanminus make libbio-samtools-perl -y --no-install-recommends
RUN cpanm HTML::Table Bio::DB::Sam
